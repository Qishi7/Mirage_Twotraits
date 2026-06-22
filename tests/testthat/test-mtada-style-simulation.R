test_that("mTADA state probabilities preserve single-trait margins", {
  pi <- mtada_state_probabilities(0.05, 0.03, 0.01)
  expect_equal(sum(pi), 1)
  expect_equal(unname(pi["10"] + pi["11"]), 0.05)
  expect_equal(unname(pi["01"] + pi["11"]), 0.03)
})

test_that("mTADA-style de novo simulation is reproducible and aligned", {
  sim1 <- simulate_mtada_style(n_genes = 40, seed = 12)
  sim2 <- simulate_mtada_style(n_genes = 40, seed = 12)
  expect_equal(sim1$truth, sim2$truth)
  expect_equal(sim1$trait1, sim2$trait1)
  expect_equal(sort(unique(sim1$trait1$Gene)), sort(sim1$truth$Gene))
  expect_true(all(is.finite(sim1$trait1$log_bf)))
  expect_true(all(sim1$trait1$No.case >= 0))
})

test_that("de novo Bayes factor matches the TADA analytic formula", {
  x <- c(0, 1, 2)
  n <- 1000
  mu <- c(1e-6, 2e-6, 3e-6)
  gamma_mean <- 12
  beta <- 0.8
  lambda <- 2 * n * mu
  expected <- log(
    dnbinom(
      x,
      size = gamma_mean * beta,
      prob = beta / (beta + lambda)
    ) / dpois(x, lambda)
  )
  expect_equal(
    log_bf_denovo_tada(x, n, mu, gamma_mean, beta),
    expected
  )
})

test_that("mTADA-style data run through the four-state MIRAGE engine", {
  sim <- simulate_mtada_style(
    n_genes = 120,
    pi_trait1 = 0.10,
    pi_trait2 = 0.08,
    pi_pleiotropic = 0.04,
    gamma_mean_trait1 = c("MiD" = 30, "LoF" = 50),
    gamma_mean_trait2 = c("MiD" = 25, "LoF" = 40),
    seed = 15
  )
  result <- fit_mtada_style_mirage(sim, max_iter = 100)
  expect_equal(
    unname(rowSums(
      result$fit$posterior[, paste0("tau_", c("00", "10", "01", "11"))]
    )),
    rep(1, nrow(result$fit$posterior)),
    tolerance = 1e-8
  )
  expect_equal(sum(result$fit$parameters$pi), 1, tolerance = 1e-8)
  expect_true(all(result$metrics$auc$auc >= 0 & result$metrics$auc$auc <= 1))
  expect_true(all(result$metrics$bfdr$observed_fdp >= 0))
  expect_true(all(result$metrics$bfdr$observed_fdp <= 1))
  expect_true(result$fit$diagnostics$loglik_non_decreasing)
})
