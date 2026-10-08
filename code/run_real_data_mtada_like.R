#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

root <- "/Users/bayeslab/Desktop/Mirage_t"
combine_dir <- file.path(root, "explorations/asc_schema_six_group_combine_in_analysis")
outdir <- file.path(combine_dir, "mtada_like_real_data")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(combine_dir, "gene_category_counts.csv")
validation_path <- file.path(
  combine_dir,
  "calibration_explore",
  "six_groups_joint_fixed_schema_delta_0.05_posterior_with_symbols_validation.csv"
)

categories <- c(
  "lof_loeuf_decile_1",
  "lof_loeuf_decile_2_3",
  "lof_loeuf_decile_4_10",
  "missense_mpc_ge_2",
  "missense_mpc_1_2",
  "missense_mpc_lt_1"
)

gamma_by_category <- c(
  lof_loeuf_decile_1 = 6,
  lof_loeuf_decile_2_3 = 6,
  lof_loeuf_decile_4_10 = 6,
  missense_mpc_ge_2 = 3,
  missense_mpc_1_2 = 3,
  missense_mpc_lt_1 = 3
)

normalize_pi <- function(pi) {
  pi <- pmax(pi, .Machine$double.xmin)
  pi / sum(pi)
}

compute_joint_quantities <- function(log_bf1, log_bf2, pi) {
  w00 <- rep(log(pi[["00"]]), length(log_bf1))
  w10 <- log(pi[["10"]]) + log_bf1
  w01 <- log(pi[["01"]]) + log_bf2
  w11 <- log(pi[["11"]]) + log_bf1 + log_bf2
  max_w <- pmax(w00, w10, w01, w11)
  denom <- exp(w00 - max_w) + exp(w10 - max_w) +
    exp(w01 - max_w) + exp(w11 - max_w)
  log_denom <- max_w + log(denom)
  tau <- cbind(
    "00" = exp(w00 - log_denom),
    "10" = exp(w10 - log_denom),
    "01" = exp(w01 - log_denom),
    "11" = exp(w11 - log_denom)
  )
  list(tau = tau, loglik = sum(log_denom))
}

update_pi_fixed_delta2 <- function(tau, fixed_delta2) {
  m <- colMeans(tau)
  out <- setNames(numeric(4), c("00", "10", "01", "11"))
  no_scz <- m[["00"]] + m[["10"]]
  scz <- m[["01"]] + m[["11"]]
  out[c("00", "10")] <- (1 - fixed_delta2) * m[c("00", "10")] / no_scz
  out[c("01", "11")] <- fixed_delta2 * m[c("01", "11")] / scz
  normalize_pi(out)
}

fit_four_state <- function(log_bf1, log_bf2, fixed_delta2 = NULL,
                           pi_init = c("00" = 0.90, "10" = 0.04,
                                       "01" = 0.04, "11" = 0.02),
                           max_iter = 500, tol = 1e-8) {
  pi <- normalize_pi(pi_init)
  trace <- vector("list", max_iter)
  for (iter in seq_len(max_iter)) {
    q <- compute_joint_quantities(log_bf1, log_bf2, pi)
    pi_new <- if (is.null(fixed_delta2)) {
      normalize_pi(colMeans(q$tau))
    } else {
      update_pi_fixed_delta2(q$tau, fixed_delta2)
    }
    trace[[iter]] <- data.table(iter = iter, loglik = q$loglik, t(pi_new))
    if (max(abs(pi_new - pi)) < tol) {
      pi <- pi_new
      break
    }
    pi <- pi_new
  }
  q <- compute_joint_quantities(log_bf1, log_bf2, pi)
  list(pi = pi, tau = q$tau, loglik = q$loglik, n_iter = iter,
       trace = rbindlist(trace[seq_len(iter)], fill = TRUE))
}

add_bfdr <- function(dt, pp_col, out_col) {
  setorderv(dt, pp_col, order = -1)
  dt[, (out_col) := cumsum(1 - get(pp_col)) / seq_len(.N)]
  invisible(dt)
}

mtada_like_log_bf <- function(d) {
  d <- copy(d)
  d[, expected_null_case := ((ctrl_ac + 0.5) / (ctrl_an_sum + 1)) * case_an_sum]
  d[, expected_alt_case := expected_null_case * gamma_by_category[category]]
  d[, log_bf_category := dpois(case_ac, pmax(expected_alt_case, .Machine$double.xmin), log = TRUE) -
      dpois(case_ac, pmax(expected_null_case, .Machine$double.xmin), log = TRUE)]
  d
}

message("Reading gene-category counts...")
counts <- fread(counts_path)
counts <- counts[category %in% categories]
counts <- mtada_like_log_bf(counts)

genes <- sort(unique(counts$gene_id))
make_gene_log_bf <- function(trait_label) {
  x <- counts[trait == trait_label]
  y <- x[, .(
    log_B = sum(log_bf_category),
    case_ac = sum(case_ac),
    ctrl_ac = sum(ctrl_ac),
    expected_null_case = sum(expected_null_case),
    variant_count = sum(variant_count)
  ), by = gene_id]
  out <- data.table(Gene = genes)
  out <- merge(out, y, by.x = "Gene", by.y = "gene_id", all.x = TRUE)
  for (col in c("log_B", "case_ac", "ctrl_ac", "expected_null_case", "variant_count")) {
    out[is.na(get(col)), (col) := 0]
  }
  out
}

asc <- make_gene_log_bf("ASC")
schema <- make_gene_log_bf("SCHEMA")

run_and_write <- function(label, fixed_delta2 = NULL) {
  fit <- fit_four_state(asc$log_B, schema$log_B, fixed_delta2 = fixed_delta2)
  posterior <- data.table(
    Gene = genes,
    tau_00 = fit$tau[, "00"],
    tau_10 = fit$tau[, "10"],
    tau_01 = fit$tau[, "01"],
    tau_11 = fit$tau[, "11"],
    PP_trait1_autism = fit$tau[, "10"] + fit$tau[, "11"],
    PP_trait2_schema = fit$tau[, "01"] + fit$tau[, "11"],
    PP_shared = fit$tau[, "11"],
    log_B_trait1_autism = asc$log_B,
    log_B_trait2_schema = schema$log_B,
    asc_case_ac = asc$case_ac,
    asc_ctrl_ac = asc$ctrl_ac,
    asc_expected_null_case = asc$expected_null_case,
    schema_case_ac = schema$case_ac,
    schema_ctrl_ac = schema$ctrl_ac,
    schema_expected_null_case = schema$expected_null_case
  )
  add_bfdr(posterior, "PP_trait1_autism", "FDR_trait1_autism")
  add_bfdr(posterior, "PP_trait2_schema", "FDR_trait2_schema")
  add_bfdr(posterior, "PP_shared", "FDR_shared")
  if (file.exists(validation_path)) {
    val <- fread(validation_path)[, .(Gene, Gene_symbol, SFARI_score1, SCHEMA_current_exomewide_12)]
    posterior <- merge(posterior, val, by = "Gene", all.x = TRUE)
  }
  fwrite(posterior, file.path(outdir, paste0(label, "_posterior.csv")))
  fwrite(fit$trace, file.path(outdir, paste0(label, "_em_trace.csv")))
  fwrite(data.table(parameter = paste0("pi_", names(fit$pi)), estimate = unname(fit$pi)),
         file.path(outdir, paste0(label, "_parameters.csv")))
  summary <- data.table(
    scenario = label,
    n_genes = nrow(posterior),
    pi00 = fit$pi[["00"]],
    pi10 = fit$pi[["10"]],
    pi01 = fit$pi[["01"]],
    pi11 = fit$pi[["11"]],
    n_iter = fit$n_iter,
    autism_fdr05 = sum(posterior$FDR_trait1_autism <= 0.05),
    schema_fdr05 = sum(posterior$FDR_trait2_schema <= 0.05),
    shared_fdr05 = sum(posterior$FDR_shared <= 0.05),
    autism_pp08 = sum(posterior$PP_trait1_autism >= 0.8),
    schema_pp08 = sum(posterior$PP_trait2_schema >= 0.8),
    shared_pp08 = sum(posterior$PP_shared >= 0.8),
    sfari_score1_in_autism_fdr05 = if ("SFARI_score1" %in% names(posterior)) {
      sum(posterior$FDR_trait1_autism <= 0.05 & posterior$SFARI_score1 == TRUE, na.rm = TRUE)
    } else NA_integer_,
    schema_known_in_schema_fdr05 = if ("SCHEMA_current_exomewide_12" %in% names(posterior)) {
      sum(posterior$FDR_trait2_schema <= 0.05 & posterior$SCHEMA_current_exomewide_12 == TRUE, na.rm = TRUE)
    } else NA_integer_
  )
  top_autism <- posterior[order(FDR_trait1_autism)][1:50]
  top_schema <- posterior[order(FDR_trait2_schema)][1:50]
  top_shared <- posterior[order(FDR_shared)][1:50]
  fwrite(top_autism, file.path(outdir, paste0(label, "_top50_autism.csv")))
  fwrite(top_schema, file.path(outdir, paste0(label, "_top50_schema.csv")))
  fwrite(top_shared, file.path(outdir, paste0(label, "_top50_shared.csv")))
  summary
}

summary <- rbindlist(list(
  run_and_write("mtada_like_free_pi"),
  run_and_write("mtada_like_fixed_schema_delta_0.05", fixed_delta2 = 0.05)
))
fwrite(summary, file.path(outdir, "mtada_like_real_data_summary.csv"))

print(summary)
