test_that("pleiotropic genes receive high posterior probability in informative simulations", {
  sim <- simulate_two_trait_mirage(
    n_genes = 200,
    variants_per_gene = 10,
    pi = c("00" = 0.65, "10" = 0.10, "01" = 0.10, "11" = 0.15),
    eta_trait1 = c("1" = 0.55, "2" = 0.15),
    eta_trait2 = c("1" = 0.55, "2" = 0.15),
    signal_mean_trait1 = 4,
    signal_mean_trait2 = 4,
    null_mean = -2.0,
    seed = 5
  )
  fit <- mirage_two_trait(
    sim$trait1, sim$trait2,
    500, 500, 500, 500,
    log_bf_col_trait1 = "log_bf",
    log_bf_col_trait2 = "log_bf",
    eta_init_trait1 = c("1" = 0.3, "2" = 0.1),
    eta_init_trait2 = c("1" = 0.3, "2" = 0.1),
    pi_init = c("00" = 0.70, "10" = 0.08, "01" = 0.08, "11" = 0.14),
    max_iter = 100
  )
  truth <- sim$truth[match(fit$posterior$Gene, sim$truth$Gene), ]
  mean_pleio <- mean(fit$posterior$PP_pleiotropy[truth$state == "11"])
  mean_non <- mean(fit$posterior$PP_pleiotropy[truth$state != "11"])
  expect_gt(mean_pleio, 0.55)
  expect_gt(mean_pleio, mean_non + 0.25)
})

test_that("eta and pi recover generating values in informative simulations", {
  true_pi <- c("00" = 0.60, "10" = 0.15, "01" = 0.10, "11" = 0.15)
  true_eta1 <- c("1" = 0.50, "2" = 0.10)
  true_eta2 <- c("1" = 0.45, "2" = 0.12)
  sim <- simulate_two_trait_mirage(
    n_genes = 300,
    variants_per_gene = 12,
    pi = true_pi,
    eta_trait1 = true_eta1,
    eta_trait2 = true_eta2,
    signal_mean_trait1 = 4,
    signal_mean_trait2 = 4,
    null_mean = -3.0,
    seed = 6
  )
  fit_pi <- mirage_two_trait(
    sim$trait1, sim$trait2,
    500, 500, 500, 500,
    log_bf_col_trait1 = "log_bf",
    log_bf_col_trait2 = "log_bf",
    fixed_eta_trait1 = true_eta1,
    fixed_eta_trait2 = true_eta2,
    pi_init = c("00" = 0.65, "10" = 0.12, "01" = 0.10, "11" = 0.13),
    max_iter = 120
  )
  fit_eta <- mirage_two_trait(
    sim$trait1, sim$trait2,
    500, 500, 500, 500,
    log_bf_col_trait1 = "log_bf",
    log_bf_col_trait2 = "log_bf",
    fixed_pi = true_pi,
    eta_init_trait1 = c("1" = 0.3, "2" = 0.05),
    eta_init_trait2 = c("1" = 0.3, "2" = 0.05),
    max_iter = 120
  )
  expect_lt(max(abs(fit_pi$parameters$pi - true_pi)), 0.18)
  expect_lt(max(abs(fit_eta$parameters$eta_trait1 - true_eta1)), 0.50)
  expect_lt(max(abs(fit_eta$parameters$eta_trait2 - true_eta2)), 0.50)
})

test_that("Bayesian FDR is calibrated in a strong simulation smoke check", {
  sim <- simulate_two_trait_mirage(
    n_genes = 250,
    variants_per_gene = 10,
    pi = c("00" = 0.70, "10" = 0.10, "01" = 0.08, "11" = 0.12),
    eta_trait1 = c("1" = 0.50, "2" = 0.12),
    eta_trait2 = c("1" = 0.50, "2" = 0.12),
    signal_mean_trait1 = 4,
    signal_mean_trait2 = 4,
    null_mean = -2.0,
    seed = 7
  )
  fit <- mirage_two_trait(
    sim$trait1, sim$trait2,
    500, 500, 500, 500,
    log_bf_col_trait1 = "log_bf",
    log_bf_col_trait2 = "log_bf",
    eta_init_trait1 = c("1" = 0.3, "2" = 0.08),
    eta_init_trait2 = c("1" = 0.3, "2" = 0.08),
    pi_init = c("00" = 0.72, "10" = 0.08, "01" = 0.08, "11" = 0.12),
    max_iter = 120
  )
  bfdr <- fit$bfdr$pleiotropy
  selected <- bfdr$Gene[bfdr$bayesian_fdr <= 0.2]
  truth <- sim$truth[match(selected, sim$truth$Gene), ]
  if (length(selected) > 0) {
    fdp <- mean(truth$state != "11")
    expect_lt(fdp, 0.45)
  }
})
