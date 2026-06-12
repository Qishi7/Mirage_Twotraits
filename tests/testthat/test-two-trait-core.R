test_that("posterior responsibilities and pi are normalized", {
  obj <- toy_fit()
  fit <- obj$fit
  tau <- as.matrix(fit$posterior[, c("tau_00", "tau_10", "tau_01", "tau_11")])
  expect_true(all(is.finite(tau)))
  expect_lt(max(abs(rowSums(tau) - 1)), 1e-10)
  expect_lt(abs(sum(fit$parameters$pi) - 1), 1e-10)
  expect_true(all(fit$parameters$pi >= 0))
})

test_that("EM observed log likelihood is nondecreasing", {
  fit <- toy_fit()$fit
  diffs <- diff(fit$traces$loglik)
  expect_true(all(diffs >= -1e-7))
  expect_true(fit$diagnostics$loglik_non_decreasing)
})

test_that("fixed independent prior reproduces separate single-trait posterior", {
  sim <- simulate_two_trait_mirage(n_genes = 60, variants_per_gene = 5, seed = 2)
  pi <- c("00" = 0.72, "10" = 0.18, "01" = 0.08, "11" = 0.02)
  # This pi factorizes with delta1 = 0.20 and delta2 = 0.10.
  fit <- mirage_two_trait(
    sim$trait1, sim$trait2,
    500, 500, 500, 500,
    log_bf_col_trait1 = "log_bf",
    log_bf_col_trait2 = "log_bf",
    fixed_pi = pi,
    fixed_eta_trait1 = c("1" = 0.20, "2" = 0.05),
    fixed_eta_trait2 = c("1" = 0.20, "2" = 0.05),
    max_iter = 3
  )
  p1 <- single_trait_posterior_from_log_bf(fit$gene_bf$log_B_trait1, 0.20)
  p2 <- single_trait_posterior_from_log_bf(fit$gene_bf$log_B_trait2, 0.10)
  expect_lt(max(abs(fit$posterior$PP_trait1 - p1)), 1e-10)
  expect_lt(max(abs(fit$posterior$PP_trait2 - p2)), 1e-10)
})

test_that("uninformative second trait does not change trait 1 under independent pi", {
  sim <- simulate_two_trait_mirage(n_genes = 60, variants_per_gene = 5, seed = 3)
  sim$trait2$log_bf <- 0
  pi <- c("00" = 0.72, "10" = 0.18, "01" = 0.08, "11" = 0.02)
  fit <- mirage_two_trait(
    sim$trait1, sim$trait2,
    500, 500, 500, 500,
    log_bf_col_trait1 = "log_bf",
    log_bf_col_trait2 = "log_bf",
    fixed_pi = pi,
    fixed_eta_trait1 = c("1" = 0.20, "2" = 0.05),
    fixed_eta_trait2 = c("1" = 0.20, "2" = 0.05),
    max_iter = 3
  )
  p1 <- single_trait_posterior_from_log_bf(fit$gene_bf$log_B_trait1, 0.20)
  expect_lt(max(abs(fit$posterior$PP_trait1 - p1)), 1e-10)
})

test_that("swapping traits swaps 10 and 01 states", {
  sim <- simulate_two_trait_mirage(n_genes = 50, variants_per_gene = 5, seed = 4)
  args <- list(
    n1_trait1 = 500, n0_trait1 = 500, n1_trait2 = 500, n0_trait2 = 500,
    log_bf_col_trait1 = "log_bf", log_bf_col_trait2 = "log_bf",
    fixed_pi = c("00" = 0.80, "10" = 0.07, "01" = 0.09, "11" = 0.04),
    fixed_eta_trait1 = c("1" = 0.20, "2" = 0.05),
    fixed_eta_trait2 = c("1" = 0.25, "2" = 0.07),
    max_iter = 3
  )
  fit12 <- do.call(mirage_two_trait, c(list(data_trait1 = sim$trait1, data_trait2 = sim$trait2), args))
  args_swapped <- args
  args_swapped$fixed_pi <- c("00" = 0.80, "10" = 0.09, "01" = 0.07, "11" = 0.04)
  args_swapped$fixed_eta_trait1 <- args$fixed_eta_trait2
  args_swapped$fixed_eta_trait2 <- args$fixed_eta_trait1
  fit21 <- do.call(mirage_two_trait, c(list(data_trait1 = sim$trait2, data_trait2 = sim$trait1), args_swapped))
  idx <- match(fit12$posterior$Gene, fit21$posterior$Gene)
  expect_lt(max(abs(fit12$posterior$tau_10 - fit21$posterior$tau_01[idx])), 1e-10)
  expect_lt(max(abs(fit12$posterior$tau_01 - fit21$posterior$tau_10[idx])), 1e-10)
  expect_lt(max(abs(fit12$posterior$tau_11 - fit21$posterior$tau_11[idx])), 1e-10)
})

test_that("invalid inputs return informative errors", {
  bad <- data.frame(Gene = "G1", No.case = -1, No.contr = 0, category = "1")
  good <- data.frame(Gene = "G1", No.case = 0, No.contr = 0, category = "1", log_bf = 0)
  expect_error(
    mirage_two_trait(bad, good, 10, 10, 10, 10, log_bf_col_trait2 = "log_bf"),
    "counts must be nonnegative"
  )
  expect_error(
    mirage_two_trait(good, good, 10, 10, 10, 10, log_bf_col_trait1 = "missing"),
    "log_bf_col"
  )
})

