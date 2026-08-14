.libPaths(c(normalizePath("Rlib", mustWork = TRUE), .libPaths()))

source("R/two_trait_mirage.R")
source("R/mtada_style_simulation.R")
source("R/mtada_full_grid.R")

output_dir <- Sys.getenv(
  "MIRAGE_SINGLE_STYLE_OUTPUT",
  "explorations/mirage_single_style_vs_mtada/results"
)
n_reps <- as.integer(Sys.getenv("MIRAGE_SINGLE_STYLE_REPS", "10"))
n_genes <- as.integer(Sys.getenv("MIRAGE_SINGLE_STYLE_GENES", "500"))
variants_per_gene <- as.integer(Sys.getenv("MIRAGE_SINGLE_STYLE_VARIANTS", "15"))
max_iter <- as.integer(Sys.getenv("MIRAGE_SINGLE_STYLE_MAX_ITER", "100"))
seed_start <- as.integer(Sys.getenv("MIRAGE_SINGLE_STYLE_SEED_START", "20260813"))
bf_mode <- match.arg(
  Sys.getenv("MIRAGE_SINGLE_STYLE_BF_MODE", "fast"),
  choices = c("fast", "integrated")
)
simulation_design <- match.arg(
  Sys.getenv("MIRAGE_SINGLE_STYLE_DESIGN", "paper_like"),
  choices = c("paper_like", "control_informative")
)
run_original <- tolower(Sys.getenv("MIRAGE_RUN_ORIGINAL_MTADA", "true")) %in%
  c("true", "1", "yes")
vb_iterations <- as.integer(Sys.getenv("MIRAGE_MTADA_VB_ITER", "3000"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "replicates"), recursive = TRUE, showWarnings = FALSE)

state_pi <- c("00" = 0.94, "10" = 0.04, "01" = 0.01, "11" = 0.01)
delta1 <- unname(state_pi[["10"]] + state_pi[["11"]])
delta2 <- unname(state_pi[["01"]] + state_pi[["11"]])
independent_pi <- c(
  "00" = (1 - delta1) * (1 - delta2),
  "10" = delta1 * (1 - delta2),
  "01" = (1 - delta1) * delta2,
  "11" = delta1 * delta2
)
eta1 <- c("C1" = 0.05, "C2" = 0.20, "C3" = 0.50)
eta2 <- c("C1" = 0.03, "C2" = 0.15, "C3" = 0.40)
rr_mean <- c("C1" = 3, "C2" = 3, "C3" = 5)
rr_sigma <- c("C1" = 1, "C2" = 1, "C3" = 1)
cat_probs <- c("C1" = 0.60, "C2" = 0.30, "C3" = 0.10)
background_prob <- as.numeric(Sys.getenv("MIRAGE_BACKGROUND_PROB", "0.18"))
background_multiplier <- as.numeric(Sys.getenv("MIRAGE_BACKGROUND_MULTIPLIER", "8"))
n_case1 <- n_control1 <- 3000L
n_case2 <- n_control2 <- 3000L

simulate_variant_level <- function(seed) {
  set.seed(seed)
  genes <- sprintf("GENE%05d", seq_len(n_genes))
  state_counts <- floor(n_genes * state_pi)
  remainder <- n_genes - sum(state_counts)
  if (remainder > 0) {
    add_order <- names(sort(n_genes * state_pi - state_counts, decreasing = TRUE))
    state_counts[add_order[seq_len(remainder)]] <-
      state_counts[add_order[seq_len(remainder)]] + 1L
  }
  states <- sample(rep(names(state_counts), state_counts))
  active1 <- states %in% c("10", "11")
  active2 <- states %in% c("01", "11")

  make_trait <- function(trait, active, eta) {
    rows <- vector("list", n_genes)
    for (i in seq_len(n_genes)) {
      cats <- sample(names(cat_probs), variants_per_gene, replace = TRUE,
                     prob = cat_probs)
      q <- pmin(rbeta(variants_per_gene, shape1 = 0.35, shape2 = 80) * 0.04,
                0.01)
      background_variant <- rep(FALSE, variants_per_gene)
      q_count <- q
      if (simulation_design == "control_informative") {
        background_variant <- rbinom(
          variants_per_gene, 1, background_prob
        ) == 1
        multiplier <- ifelse(background_variant, background_multiplier, 1)
        q_count <- pmin(q * multiplier, 0.05)
      }
      z <- if (active[i]) {
        rbinom(variants_per_gene, 1, eta[cats])
      } else {
        rep(0L, variants_per_gene)
      }
      rr <- rep(1, variants_per_gene)
      idx <- z == 1
      if (any(idx)) {
        rr[idx] <- rgamma(
          sum(idx),
          shape = rr_mean[cats[idx]] * rr_sigma[cats[idx]],
          rate = rr_sigma[cats[idx]]
        )
      }
      n_case <- if (trait == 1) n_case1 else n_case2
      n_control <- if (trait == 1) n_control1 else n_control2
      case_count <- rpois(variants_per_gene, 2 * n_case * q_count * rr)
      control_count <- rpois(variants_per_gene, 2 * n_control * q_count)
      total_count <- case_count + control_count
      if (bf_mode == "fast") {
        p_null <- n_case / (n_case + n_control)
        p_alt <- rr_mean[cats] * n_case / (rr_mean[cats] * n_case + n_control)
        log_bf <- dbinom(case_count, total_count, p_alt, log = TRUE) -
          dbinom(case_count, total_count, p_null, log = TRUE)
      } else {
        log_bf <- NA_real_
      }
      rows[[i]] <- data.frame(
        ID = paste0("t", trait, "_", genes[i], "_v", seq_len(variants_per_gene)),
        Gene = genes[i],
        No.case = case_count,
        No.contr = control_count,
        category = cats,
        log_bf = log_bf,
        mutation_rate = q,
        count_rate = q_count,
        background_variant = background_variant,
        true_z = z,
        true_rr = rr,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }

  list(
    trait1 = make_trait(1, active1, eta1),
    trait2 = make_trait(2, active2, eta2),
    truth = data.frame(
      Gene = genes,
      state = states,
      trait1_risk = active1,
      trait2_risk = active2,
      pleiotropic = states == "11",
      stringsAsFactors = FALSE
    ),
    parameters = list(
      pi = state_pi,
      eta_trait1 = eta1,
      eta_trait2 = eta2,
      rr_mean = rr_mean,
      rr_sigma = rr_sigma,
      simulation_design = simulation_design,
      background_prob = background_prob,
      background_multiplier = background_multiplier
    )
  )
}

posterior_from_fit <- function(fit) {
  fit$posterior[, c("Gene", "PP_trait1", "PP_trait2", "PP_pleiotropy")]
}

aggregate_for_mtada <- function(trait_data) {
  counts <- aggregate(No.case ~ Gene + category, trait_data, sum)
  rates <- aggregate(mutation_rate ~ Gene + category, trait_data, sum)
  genes <- sort(unique(trait_data$Gene))
  cats <- names(cat_probs)
  count_mat <- matrix(0L, nrow = length(genes), ncol = length(cats),
                      dimnames = list(genes, cats))
  rate_mat <- matrix(.Machine$double.xmin, nrow = length(genes),
                     ncol = length(cats), dimnames = list(genes, cats))
  count_mat[cbind(counts$Gene, counts$category)] <- counts$No.case
  rate_mat[cbind(rates$Gene, rates$category)] <- pmax(
    rates$mutation_rate, .Machine$double.xmin
  )
  list(counts = count_mat, rates = rate_mat)
}

fit_mtada_aggregated <- function(sim, compiled, seed) {
  env <- compiled$env
  agg1 <- aggregate_for_mtada(sim$trait1)
  agg2 <- aggregate_for_mtada(sim$trait2)
  genes <- rownames(agg1$counts)
  cats <- colnames(agg1$counts)
  stan_data <- list(
    NN = length(genes),
    NCdn1 = length(cats),
    Ndn1 = rep(n_case1, length(cats)),
    hyperGammaMeanDN1 = unname(rr_mean[cats]),
    NCdn2 = length(cats),
    Ndn2 = rep(n_case2, length(cats)),
    hyperGammaMeanDN2 = unname(rr_mean[cats]),
    dataDN1 = agg1$counts,
    mutRate1 = agg1$rates,
    dataDN2 = agg2$counts,
    mutRate2 = agg2$rates,
    betaPars = c(6.7771073, -1.7950864, -0.2168248),
    lowerGamma = 1,
    lowerBeta = 1,
    hyperBetaDN01 = rep(1, length(cats)),
    hyperBetaDN02 = rep(1, length(cats)),
    adjustHyperBeta = 0L,
    pi01 = delta1,
    pi02 = delta2
  )
  vb_fit <- rstan::vb(
    object = compiled$model,
    data = stan_data,
    pars = c("p12", "gammaMeanDN1"),
    init = list(p12 = min(state_pi[["11"]], 0.9 * min(delta1, delta2))),
    iter = vb_iterations,
    seed = seed,
    refresh = 0
  )
  estimates <- tryCatch(
    env$estimatePars(c("p12", "gammaMeanDN1[1]"), vb_fit),
    error = function(e) e
  )
  if (inherits(estimates, "error")) {
    draws <- as.data.frame(vb_fit)[["p12"]]
    density_fit <- density(draws, from = 0, to = min(delta1, delta2), cut = 0)
    pi11 <- density_fit$x[which.max(density_fit$y)]
    estimator <- paste0("density_fallback: ", conditionMessage(estimates))
  } else {
    pi11 <- as.numeric(estimates["p12", "Mode"])
    estimator <- "mTADA_locfit_mode"
  }
  if (pi11 < 1e-4) pi11 <- 0
  priors <- c(1 - delta1 - delta2 + pi11, pi11, delta1 - pi11, delta2 - pi11)
  gamma_states <- rbind(
    rep(1, 2 * length(cats)),
    c(rr_mean[cats], rr_mean[cats]),
    c(rr_mean[cats], rep(1, length(cats))),
    c(rep(1, length(cats)), rr_mean[cats])
  )
  pp <- env$posProb.dn(
    dnData = data.frame(agg1$counts, agg2$counts),
    muAll = data.frame(agg1$rates, agg2$rates),
    gamma.mean.dn = gamma_states,
    Ndn = c(stan_data$Ndn1, stan_data$Ndn2),
    prob0 = priors,
    beta.dn = matrix(1, nrow = 4, ncol = 2 * length(cats))
  )$PP
  out <- data.frame(
    Gene = genes,
    PP_trait1 = pp[, "FIRST"] + pp[, "BOTH"],
    PP_trait2 = pp[, "SECOND"] + pp[, "BOTH"],
    PP_pleiotropy = pp[, "BOTH"],
    stringsAsFactors = FALSE
  )
  attr(out, "pi11_estimate") <- pi11
  attr(out, "estimator") <- estimator
  out
}

parameter_rows <- function(seed, replicate, method, pi_est, eta1_est = NULL,
                           eta2_est = NULL) {
  pi_truth <- c("pi00" = state_pi[["00"]], "pi10" = state_pi[["10"]],
                "pi01" = state_pi[["01"]], "pi11" = state_pi[["11"]])
  pi_out <- data.frame(
    seed = seed,
    replicate = replicate,
    method = method,
    parameter = names(pi_truth),
    truth = unname(pi_truth),
    estimate = unname(pi_est[names(pi_truth)]),
    stringsAsFactors = FALSE
  )
  rows <- list(pi_out)
  if (!is.null(eta1_est)) {
    rows[[length(rows) + 1]] <- data.frame(
      seed = seed, replicate = replicate, method = method,
      parameter = paste0("eta1_", names(eta1)),
      truth = unname(eta1),
      estimate = unname(eta1_est[names(eta1)]),
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(eta2_est)) {
    rows[[length(rows) + 1]] <- data.frame(
      seed = seed, replicate = replicate, method = method,
      parameter = paste0("eta2_", names(eta2)),
      truth = unname(eta2),
      estimate = unname(eta2_est[names(eta2)]),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out$abs_error <- abs(out$estimate - out$truth)
  out
}

collect_results <- function() {
  paths <- list.files(file.path(output_dir, "replicates"),
                      pattern = "^rep_[0-9]+\\.rds$", full.names = TRUE)
  if (!length(paths)) {
    return(invisible(list(metrics = data.frame(), parameters = data.frame())))
  }
  results <- lapply(paths, readRDS)
  metrics <- do.call(rbind, lapply(results, function(x) {
    cbind(x$meta[rep(1, nrow(x$metrics)), , drop = FALSE], x$metrics)
  }))
  params <- do.call(rbind, lapply(results, `[[`, "parameters"))
  write.csv(metrics, file.path(output_dir, "all_metrics.csv"), row.names = FALSE)
  write.csv(params, file.path(output_dir, "parameter_learning.csv"), row.names = FALSE)
  auc <- metrics[metrics$metric == "auc", , drop = FALSE]
  auc_summary <- do.call(rbind, lapply(
    split(auc, paste(auc$method, auc$target, sep = "|")),
    function(z) data.frame(method = z$method[1], target = z$target[1],
                           mean_auc = mean(z$value), sd_auc = sd(z$value),
                           min_auc = min(z$value), max_auc = max(z$value),
                           n = nrow(z), stringsAsFactors = FALSE)
  ))
  rownames(auc_summary) <- NULL
  write.csv(auc_summary, file.path(output_dir, "auc_summary.csv"), row.names = FALSE)
  param_summary <- do.call(rbind, lapply(
    split(params, paste(params$method, params$parameter, sep = "|")),
    function(z) data.frame(method = z$method[1], parameter = z$parameter[1],
                           mean_estimate = mean(z$estimate, na.rm = TRUE),
                           sd_estimate = sd(z$estimate, na.rm = TRUE),
                           mean_abs_error = mean(z$abs_error, na.rm = TRUE),
                           sd_abs_error = sd(z$abs_error, na.rm = TRUE),
                           n = nrow(z), stringsAsFactors = FALSE)
  ))
  rownames(param_summary) <- NULL
  write.csv(param_summary, file.path(output_dir, "parameter_summary.csv"),
            row.names = FALSE)
  invisible(list(metrics = metrics, parameters = params))
}

compiled <- NULL
if (run_original) {
  compiled <- compile_original_mtada(".tmp_mtada_reference/script/mTADA.R")
  write.csv(data.frame(method = "mtada_aggregated",
                       compile_seconds = compiled$compile_seconds),
            file.path(output_dir, "compile_runtime.csv"), row.names = FALSE)
}

for (replicate in seq_len(n_reps)) {
  seed <- seed_start + replicate - 1L
  path <- file.path(output_dir, "replicates", sprintf("rep_%03d.rds", replicate))
  if (file.exists(path)) {
    message("Skipping existing replicate ", replicate)
    next
  }
  message("Running replicate ", replicate, "/", n_reps, " seed ", seed)
  sim <- simulate_variant_level(seed)
  truth <- sim$truth
  log_bf_col <- if (bf_mode == "fast") "log_bf" else NULL

  joint <- time_call(mirage_two_trait(
    sim$trait1, sim$trait2,
    n1_trait1 = n_case1, n0_trait1 = n_control1,
    n1_trait2 = n_case2, n0_trait2 = n_control2,
    gamma_trait1 = rr_mean, sigma_trait1 = rr_sigma,
    gamma_trait2 = rr_mean, sigma_trait2 = rr_sigma,
    log_bf_col_trait1 = log_bf_col,
    log_bf_col_trait2 = log_bf_col,
    pi_init = state_pi,
    eta_init_trait1 = eta1,
    eta_init_trait2 = eta2,
    max_iter = max_iter
  ))
  joint_post <- posterior_from_fit(joint$value)
  metrics <- evaluate_posterior_method(
    "mirage_joint", joint_post, truth, joint$seconds,
    joint$value$parameters$pi[["11"]]
  )
  joint_pi <- c("pi00" = joint$value$parameters$pi[["00"]],
                "pi10" = joint$value$parameters$pi[["10"]],
                "pi01" = joint$value$parameters$pi[["01"]],
                "pi11" = joint$value$parameters$pi[["11"]])
  parameters <- parameter_rows(
    seed, replicate, "mirage_joint", joint_pi,
    joint$value$parameters$eta_trait1, joint$value$parameters$eta_trait2
  )

  sep <- time_call(mirage_two_trait(
    sim$trait1, sim$trait2,
    n1_trait1 = n_case1, n0_trait1 = n_control1,
    n1_trait2 = n_case2, n0_trait2 = n_control2,
    gamma_trait1 = rr_mean, sigma_trait1 = rr_sigma,
    gamma_trait2 = rr_mean, sigma_trait2 = rr_sigma,
    log_bf_col_trait1 = log_bf_col,
    log_bf_col_trait2 = log_bf_col,
    fixed_pi = independent_pi,
    eta_init_trait1 = eta1,
    eta_init_trait2 = eta2,
    max_iter = max_iter
  ))
  sep_post <- posterior_from_fit(sep$value)
  metrics <- rbind(metrics, evaluate_posterior_method(
    "separate_mirage", sep_post, truth, sep$seconds,
    independent_pi[["11"]]
  ))
  sep_pi <- c("pi00" = independent_pi[["00"]], "pi10" = independent_pi[["10"]],
              "pi01" = independent_pi[["01"]], "pi11" = independent_pi[["11"]])
  parameters <- rbind(parameters, parameter_rows(
    seed, replicate, "separate_mirage", sep_pi,
    sep$value$parameters$eta_trait1, sep$value$parameters$eta_trait2
  ))

  mtada_error <- ""
  if (run_original) {
    mtada <- tryCatch(
      time_call(fit_mtada_aggregated(sim, compiled, seed)),
      error = function(e) e
    )
    if (inherits(mtada, "error")) {
      mtada_error <- conditionMessage(mtada)
    } else {
      pi11 <- attr(mtada$value, "pi11_estimate")
      metrics <- rbind(metrics, evaluate_posterior_method(
        "mtada_aggregated", mtada$value, truth, mtada$seconds, pi11
      ))
      mtada_pi <- c("pi00" = 1 - delta1 - delta2 + pi11,
                    "pi10" = delta1 - pi11,
                    "pi01" = delta2 - pi11,
                    "pi11" = pi11)
      parameters <- rbind(parameters, parameter_rows(
        seed, replicate, "mtada_aggregated", mtada_pi
      ))
    }
  }

  meta <- data.frame(
    seed = seed,
    replicate = replicate,
    n_genes = n_genes,
    variants_per_gene = variants_per_gene,
    n_case1 = n_case1,
    n_control1 = n_control1,
    n_case2 = n_case2,
    n_control2 = n_control2,
    bf_mode = bf_mode,
    simulation_design = simulation_design,
    background_prob = background_prob,
    background_multiplier = background_multiplier,
    pi11_truth = state_pi[["11"]],
    mtada_error = mtada_error,
    stringsAsFactors = FALSE
  )
  saveRDS(list(meta = meta, metrics = metrics, parameters = parameters),
          paste0(path, ".tmp"))
  file.rename(paste0(path, ".tmp"), path)
  collect_results()
}

collect_results()
message("Finished exploration: ", output_dir)
