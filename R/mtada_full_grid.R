# Full mTADA paper simulation grid and method comparison.

mtada_paper_grid <- function(n_replicates = 100L) {
  rr <- list(
    high_low = list(trait1 = c("MiD" = 105, "LoF" = 29),
                    trait2 = c("MiD" = 12, "LoF" = 2)),
    high_high = list(trait1 = c("MiD" = 105, "LoF" = 29),
                     trait2 = c("MiD" = 87, "LoF" = 23)),
    low_low = list(trait1 = c("MiD" = 24, "LoF" = 6),
                   trait2 = c("MiD" = 12, "LoF" = 2)),
    low_high = list(trait1 = c("MiD" = 24, "LoF" = 6),
                    trait2 = c("MiD" = 87, "LoF" = 23))
  )
  design <- expand.grid(
    pi_pleiotropic = c(0, 0.02),
    rr_setting = names(rr),
    n_trios_trait1 = c(1000L, 2000L, 5000L),
    n_trios_trait2 = c(1000L, 2000L, 5000L),
    replicate = seq_len(as.integer(n_replicates)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  design$gamma_mean_trait1 <- lapply(
    design$rr_setting, function(x) rr[[x]]$trait1
  )
  design$gamma_mean_trait2 <- lapply(
    design$rr_setting, function(x) rr[[x]]$trait2
  )
  design$task_id <- seq_len(nrow(design))
  design$seed <- 20200610L + design$task_id
  design[, c(
    "task_id", "seed", "replicate", "pi_pleiotropic", "rr_setting",
    "n_trios_trait1", "n_trios_trait2",
    "gamma_mean_trait1", "gamma_mean_trait2"
  )]
}

load_mtada_paper_mutation_rates <- function(path) {
  if (!file.exists(path)) {
    stop("mTADA mutation-rate file was not found: ", path, call. = FALSE)
  }
  header <- names(read.table(path, header = TRUE, nrows = 1,
                             stringsAsFactors = FALSE))
  required <- c("Gene", "mut_damaging", "mut_lof")
  if (!all(required %in% header)) {
    stop("mTADA mutation-rate file is missing required columns.",
         call. = FALSE)
  }
  data <- read.table(
    path,
    header = TRUE,
    colClasses = ifelse(header %in% required, NA, "NULL"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rates <- as.matrix(data[, c("mut_damaging", "mut_lof")])
  storage.mode(rates) <- "double"
  colnames(rates) <- c("MiD", "LoF")
  rownames(rates) <- data$Gene
  if (nrow(rates) != 19358L) {
    warning("The paper data contain 19,358 genes; loaded ", nrow(rates), ".")
  }
  rates[rates <= 0] <- .Machine$double.xmin
  rates
}

time_call <- function(expr) {
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, seconds = proc.time()[["elapsed"]] - start)
}

fit_four_state_log_bf <- function(log_bf_trait1, log_bf_trait2,
                                  pi_init = c(
                                    "00" = 0.92, "10" = 0.04,
                                    "01" = 0.02, "11" = 0.02
                                  ),
                                  tolerance = 1e-6,
                                  max_iter = 500L) {
  fast_quantities <- function(pi) {
    w00 <- rep(log(pi[["00"]]), length(log_bf_trait1))
    w10 <- log(pi[["10"]]) + log_bf_trait1
    w01 <- log(pi[["01"]]) + log_bf_trait2
    w11 <- log(pi[["11"]]) + log_bf_trait1 + log_bf_trait2
    max_w <- pmax(w00, w10, w01, w11)
    denominator <- exp(w00 - max_w) + exp(w10 - max_w) +
      exp(w01 - max_w) + exp(w11 - max_w)
    log_denominator <- max_w + log(denominator)
    tau <- cbind(
      "00" = exp(w00 - log_denominator),
      "10" = exp(w10 - log_denominator),
      "01" = exp(w01 - log_denominator),
      "11" = exp(w11 - log_denominator)
    )
    list(tau = tau, loglik = sum(log_denominator))
  }
  pi <- normalize_pi(pi_init)
  trace <- numeric(max_iter)
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    quantities <- fast_quantities(pi)
    pi_new <- normalize_pi(colMeans(quantities$tau))
    trace[iter] <- quantities$loglik
    if (max(abs(pi_new - pi)) < tolerance) {
      pi <- pi_new
      converged <- TRUE
      break
    }
    pi <- pi_new
  }
  quantities <- fast_quantities(pi)
  list(
    pi = pi,
    tau = quantities$tau,
    loglik = quantities$loglik,
    converged = converged,
    n_iter = iter,
    trace_loglik = trace[seq_len(iter)]
  )
}

fit_joint_mtada_style <- function(sim) {
  genes <- sim$truth$Gene
  log_bf1 <- tapply(sim$trait1$log_bf, sim$trait1$Gene, sum)
  log_bf2 <- tapply(sim$trait2$log_bf, sim$trait2$Gene, sum)
  fit <- fit_four_state_log_bf(
    unname(log_bf1[genes]), unname(log_bf2[genes])
  )
  posterior <- data.frame(
    Gene = genes,
    PP_trait1 = fit$tau[, "10"] + fit$tau[, "11"],
    PP_trait2 = fit$tau[, "01"] + fit$tau[, "11"],
    PP_pleiotropy = fit$tau[, "11"],
    stringsAsFactors = FALSE
  )
  attr(posterior, "pi") <- fit$pi
  attr(posterior, "diagnostics") <- fit[
    c("converged", "n_iter", "loglik", "trace_loglik")
  ]
  posterior
}

evaluate_posterior_method <- function(method, posterior, truth,
                                      runtime_seconds,
                                      pi11_estimate = NA_real_,
                                      pp_threshold = 0.8,
                                      bfdr_targets = c(0.01, 0.05, 0.10, 0.20)) {
  targets <- list(
    trait1 = list(score = posterior$PP_trait1, truth = truth$trait1_risk),
    trait2 = list(score = posterior$PP_trait2, truth = truth$trait2_risk),
    pleiotropy = list(
      score = posterior$PP_pleiotropy, truth = truth$pleiotropic
    )
  )
  auc <- do.call(rbind, lapply(names(targets), function(target) {
    data.frame(
      method = method,
      target = target,
      metric = "auc",
      threshold = NA_real_,
      value = binary_auc(
        targets[[target]]$score, targets[[target]]$truth
      ),
      stringsAsFactors = FALSE
    )
  }))
  pp <- do.call(rbind, lapply(names(targets), function(target) {
    x <- selection_metrics(
      targets[[target]]$score,
      targets[[target]]$truth,
      pp_threshold
    )
    rbind(
      data.frame(method = method, target = target, metric = "pp08_power",
                 threshold = pp_threshold, value = x$power),
      data.frame(method = method, target = target, metric = "pp08_fdp",
                 threshold = pp_threshold, value = x$fdp),
      data.frame(method = method, target = target,
                 metric = "pp08_discoveries",
                 threshold = pp_threshold, value = x$discoveries)
    )
  }))
  calibration <- do.call(rbind, lapply(names(targets), function(target) {
    do.call(rbind, lapply(bfdr_targets, function(alpha) {
      x <- bfdr_metrics(
        targets[[target]]$score, targets[[target]]$truth, alpha
      )
      rbind(
        data.frame(method = method, target = target,
                   metric = "bfdr_observed_fdp",
                   threshold = alpha, value = x$observed_fdp),
        data.frame(method = method, target = target,
                   metric = "bfdr_estimate",
                   threshold = alpha, value = x$estimated_bfdr),
        data.frame(method = method, target = target,
                   metric = "bfdr_power",
                   threshold = alpha, value = x$power),
        data.frame(method = method, target = target,
                   metric = "bfdr_discoveries",
                   threshold = alpha, value = x$discoveries)
      )
    }))
  }))
  extra <- rbind(
    data.frame(method = method, target = "global", metric = "runtime_seconds",
               threshold = NA_real_, value = runtime_seconds),
    data.frame(method = method, target = "pleiotropy",
               metric = "pi11_estimate",
               threshold = NA_real_, value = pi11_estimate)
  )
  rbind(auc, pp, calibration, extra)
}

fit_separate_mtada_style <- function(sim) {
  log_bf1 <- tapply(sim$trait1$log_bf, sim$trait1$Gene, sum)
  log_bf2 <- tapply(sim$trait2$log_bf, sim$trait2$Gene, sum)
  genes <- sim$truth$Gene
  data.frame(
    Gene = genes,
    PP_trait1 = single_trait_posterior_from_log_bf(
      unname(log_bf1[genes]), sim$parameters$pi_trait1
    ),
    PP_trait2 = single_trait_posterior_from_log_bf(
      unname(log_bf2[genes]), sim$parameters$pi_trait2
    ),
    PP_pleiotropy = single_trait_posterior_from_log_bf(
      unname(log_bf1[genes]), sim$parameters$pi_trait1
    ) * single_trait_posterior_from_log_bf(
      unname(log_bf2[genes]), sim$parameters$pi_trait2
    ),
    stringsAsFactors = FALSE
  )
}

compile_original_mtada <- function(source_file) {
  if (!requireNamespace("rstan", quietly = TRUE) ||
      !requireNamespace("locfit", quietly = TRUE)) {
    stop("Original mTADA requires rstan and locfit.", call. = FALSE)
  }
  env <- new.env(parent = globalenv())
  sys.source(source_file, envir = env)
  compiled <- time_call(rstan::stan_model(model_code = env$DN2traits))
  list(env = env, model = compiled$value,
       compile_seconds = compiled$seconds)
}

fit_original_mtada <- function(sim, compiled, vb_iterations = 5000L,
                               seed = 1L) {
  env <- compiled$env
  categories <- names(sim$parameters$gamma_mean_trait1)
  genes <- sim$truth$Gene
  counts1 <- reshape(
    sim$trait1[, c("Gene", "category", "No.case")],
    idvar = "Gene", timevar = "category", direction = "wide"
  )
  counts2 <- reshape(
    sim$trait2[, c("Gene", "category", "No.case")],
    idvar = "Gene", timevar = "category", direction = "wide"
  )
  rates1 <- reshape(
    sim$trait1[, c("Gene", "category", "mutation_rate")],
    idvar = "Gene", timevar = "category", direction = "wide"
  )
  rates2 <- reshape(
    sim$trait2[, c("Gene", "category", "mutation_rate")],
    idvar = "Gene", timevar = "category", direction = "wide"
  )
  reorder_matrix <- function(data, prefix) {
    data <- data[match(genes, data$Gene), , drop = FALSE]
    out <- as.matrix(data[, paste0(prefix, ".", categories), drop = FALSE])
    storage.mode(out) <- "double"
    out
  }
  data1 <- reorder_matrix(counts1, "No.case")
  data2 <- reorder_matrix(counts2, "No.case")
  mu1 <- reorder_matrix(rates1, "mutation_rate")
  mu2 <- reorder_matrix(rates2, "mutation_rate")
  p1 <- sim$parameters$pi_trait1
  p2 <- sim$parameters$pi_trait2
  init_p12 <- max(min(sim$parameters$pi_pleiotropic, 0.9 * min(p1, p2)),
                  0.001 * min(p1, p2))
  stan_data <- list(
    NN = length(genes),
    NCdn1 = length(categories),
    Ndn1 = rep(sim$parameters$n_trios_trait1, length(categories)),
    hyperGammaMeanDN1 = unname(sim$parameters$gamma_mean_trait1),
    NCdn2 = length(categories),
    Ndn2 = rep(sim$parameters$n_trios_trait2, length(categories)),
    hyperGammaMeanDN2 = unname(sim$parameters$gamma_mean_trait2),
    dataDN1 = data1,
    mutRate1 = mu1,
    dataDN2 = data2,
    mutRate2 = mu2,
    betaPars = c(6.7771073, -1.7950864, -0.2168248),
    lowerGamma = 1,
    lowerBeta = 1,
    hyperBetaDN01 = rep(1, length(categories)),
    hyperBetaDN02 = rep(1, length(categories)),
    adjustHyperBeta = 0L,
    pi01 = p1,
    pi02 = p2
  )
  vb_fit <- rstan::vb(
    object = compiled$model,
    data = stan_data,
    pars = c("p12", "gammaMeanDN1"),
    init = list(p12 = init_p12),
    iter = as.integer(vb_iterations),
    seed = as.integer(seed),
    refresh = 0
  )
  estimator <- "mTADA_locfit_mode"
  estimates <- tryCatch(
    env$estimatePars(
      pars = c("p12", "gammaMeanDN1[1]"), mcmcResult = vb_fit
    ),
    error = function(e) e
  )
  if (inherits(estimates, "error")) {
    draws <- as.data.frame(vb_fit)[["p12"]]
    density_fit <- density(
      draws,
      from = 0,
      to = min(p1, p2),
      cut = 0
    )
    pi11 <- density_fit$x[which.max(density_fit$y)]
    estimator <- paste0(
      "bounded_density_mode_fallback: ", conditionMessage(estimates)
    )
  } else {
    pi11 <- as.numeric(estimates["p12", "Mode"])
  }
  if (pi11 < 1e-4) pi11 <- 0
  priors <- c(
    1 - p1 - p2 + pi11,
    pi11,
    p1 - pi11,
    p2 - pi11
  )
  gamma_states <- rbind(
    rep(1, 2 * length(categories)),
    c(sim$parameters$gamma_mean_trait1,
      sim$parameters$gamma_mean_trait2),
    c(sim$parameters$gamma_mean_trait1,
      rep(1, length(categories))),
    c(rep(1, length(categories)),
      sim$parameters$gamma_mean_trait2)
  )
  posterior <- env$posProb.dn(
    dnData = data.frame(data1, data2),
    muAll = data.frame(mu1, mu2),
    gamma.mean.dn = gamma_states,
    Ndn = c(stan_data$Ndn1, stan_data$Ndn2),
    prob0 = priors,
    beta.dn = matrix(1, nrow = 4, ncol = 2 * length(categories))
  )$PP
  data.frame(
    Gene = genes,
    PP_trait1 = posterior[, "FIRST"] + posterior[, "BOTH"],
    PP_trait2 = posterior[, "SECOND"] + posterior[, "BOTH"],
    PP_pleiotropy = posterior[, "BOTH"],
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    structure(pi11_estimate = pi11, pi11_estimator = estimator)
}

run_mtada_grid_task <- function(task, mutation_rates,
                                original_mtada = NULL,
                                run_original_mtada = TRUE,
                                vb_iterations = 5000L) {
  sim <- simulate_mtada_style(
    n_genes = nrow(mutation_rates),
    mutation_rates = mutation_rates,
    pi_trait1 = 0.05,
    pi_trait2 = 0.03,
    pi_pleiotropic = task$pi_pleiotropic,
    n_trios_trait1 = task$n_trios_trait1,
    n_trios_trait2 = task$n_trios_trait2,
    gamma_mean_trait1 = task$gamma_mean_trait1[[1]],
    gamma_mean_trait2 = task$gamma_mean_trait2[[1]],
    beta_trait1 = c("MiD" = 1, "LoF" = 1),
    beta_trait2 = c("MiD" = 1, "LoF" = 1),
    seed = task$seed
  )
  truth <- sim$truth

  joint <- time_call(fit_joint_mtada_style(sim))
  joint_pi <- attr(joint$value, "pi")
  joint_diagnostics <- attr(joint$value, "diagnostics")
  metrics <- evaluate_posterior_method(
    "mirage_joint", joint$value, truth, joint$seconds,
    joint_pi[["11"]]
  )

  separate <- time_call(fit_separate_mtada_style(sim))
  metrics <- rbind(metrics, evaluate_posterior_method(
    "separate", separate$value, truth, separate$seconds
  ))

  errors <- character()
  original_estimator <- NA_character_
  if (run_original_mtada && !is.null(original_mtada)) {
    original <- tryCatch(
      time_call(fit_original_mtada(
        sim, original_mtada, vb_iterations, task$seed
      )),
      error = function(e) e
    )
    if (inherits(original, "error")) {
      errors <- conditionMessage(original)
    } else {
      metrics <- rbind(metrics, evaluate_posterior_method(
        "mtada_original", original$value, truth, original$seconds,
        attr(original$value, "pi11_estimate")
      ))
      original_estimator <- attr(original$value, "pi11_estimator")
    }
  }
  meta <- task[, setdiff(names(task),
                         c("gamma_mean_trait1", "gamma_mean_trait2")),
               drop = FALSE]
  meta$pi11_truth <- task$pi_pleiotropic
  meta$joint_converged <- joint_diagnostics$converged
  meta$joint_n_iter <- joint_diagnostics$n_iter
  meta$joint_loglik_monotone <- all(
    diff(joint_diagnostics$trace_loglik) >= -1e-8
  )
  meta$original_mtada_error <- paste(errors, collapse = " | ")
  meta$original_mtada_pi11_estimator <- original_estimator
  list(meta = meta, metrics = metrics)
}

run_mtada_full_grid <- function(
    mutation_rate_file,
    output_dir = "simulation_results/mtada_full_grid",
    task_ids = NULL,
    n_replicates = 100L,
    mtada_source_file = NULL,
    run_original_mtada = TRUE,
    vb_iterations = 5000L) {
  design <- mtada_paper_grid(n_replicates)
  if (!is.null(task_ids)) {
    design <- design[design$task_id %in% task_ids, , drop = FALSE]
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  mutation_rates <- load_mtada_paper_mutation_rates(mutation_rate_file)

  original <- NULL
  compile_info <- NULL
  if (run_original_mtada) {
    if (is.null(mtada_source_file) || !file.exists(mtada_source_file)) {
      warning("Original mTADA source was not found; continuing without it.")
      run_original_mtada <- FALSE
    } else {
      original <- compile_original_mtada(mtada_source_file)
      compile_info <- data.frame(
        method = "mtada_original",
        compile_seconds = original$compile_seconds
      )
      write.csv(compile_info, file.path(output_dir, "compile_runtime.csv"),
                row.names = FALSE)
    }
  }

  for (i in seq_len(nrow(design))) {
    task <- design[i, , drop = FALSE]
    path <- file.path(
      output_dir, sprintf("task_%05d.rds", task$task_id)
    )
    if (file.exists(path)) next
    result <- tryCatch(
      run_mtada_grid_task(
        task, mutation_rates, original, run_original_mtada, vb_iterations
      ),
      error = function(e) list(
        meta = transform(
          task[, setdiff(names(task),
                         c("gamma_mean_trait1", "gamma_mean_trait2")),
               drop = FALSE],
          fatal_error = conditionMessage(e)
        ),
        metrics = data.frame()
      )
    )
    tmp <- paste0(path, ".tmp")
    saveRDS(result, tmp)
    if (!file.rename(tmp, path)) {
      stop("Could not finalize checkpoint: ", path, call. = FALSE)
    }
  }
  invisible(design)
}

collect_mtada_grid_results <- function(output_dir) {
  paths <- list.files(output_dir, pattern = "^task_[0-9]+\\.rds$",
                      full.names = TRUE)
  results <- lapply(paths, readRDS)
  metrics <- do.call(rbind, lapply(results, function(x) {
    if (!nrow(x$metrics)) return(NULL)
    cbind(x$meta[rep(1, nrow(x$metrics)), , drop = FALSE], x$metrics)
  }))
  errors <- do.call(rbind, lapply(results, function(x) x$meta))
  list(metrics = metrics, task_status = errors)
}
