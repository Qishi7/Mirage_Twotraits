convergence_data <- function(lb = c(-1, 2, -2, 3)) {
  data.frame(Gene = paste0("g", seq_along(lb)), ID = paste0("v", seq_along(lb)),
             No.case = 0, No.contr = 0, category = "c", log_bf = lb)
}
convergence_fit <- function(...) {
  d <- convergence_data()
  mirage_two_trait(d, d, 100, 100, 100, 100,
    log_bf_col_trait1 = "log_bf", log_bf_col_trait2 = "log_bf", ...)
}

test_that("iteration limit and actual stopping metrics are reported", {
  fit <- convergence_fit(max_iter = 1)
  d <- fit$diagnostics$convergence
  expect_false(d$converged)
  expect_identical(d$stop_reason, "iteration_limit")
  expect_equal(d$final_parameter_change, tail(fit$traces$max_parameter_change, 1))
  expect_equal(d$final_relative_loglik_change,
    abs(fit$traces$loglik[1] - fit$diagnostics$initial_loglik) /
      (1 + abs(fit$diagnostics$initial_loglik)))
  expect_match(paste(d$warnings, collapse = " "), "not met")
  expect_output(print(fit), "iteration_limit")
  expect_equal(fit$traces$objective, fit$traces$loglik)
  expect_equal(mirage_convergence(fit), d)
})

test_that("small EM steps at an absorbing boundary do not imply stationarity", {
  a <- convergence_data(rep(3, 4)); b <- convergence_data(rep(0, 4))
  fit <- mirage_two_trait(a, b, 100, 100, 100, 100,
    log_bf_col_trait1 = "log_bf", log_bf_col_trait2 = "log_bf",
    fixed_pi = c(.5, .5, 0, 0), eta_init_trait1 = 0, eta_init_trait2 = .2)
  d <- mirage_convergence(fit, hessian = TRUE)
  expect_true(d$converged)
  expect_false(d$stationary)
  expect_gt(d$eta_kkt_residual, 1)
  expect_match(paste(d$warnings, collapse = " "), "stationarity check failed")
  expect_identical(d$information$status, "boundary_not_regular")
  expect_true(any(d$boundary$near_boundary & d$boundary$estimated))
})

test_that("fixed parameters are excluded from stationarity and curvature checks", {
  fit <- convergence_fit(fixed_pi = c(1, 0, 0, 0),
    fixed_eta_trait1 = 0, fixed_eta_trait2 = 1)
  d <- mirage_convergence(fit, hessian = TRUE)
  expect_true(d$converged)
  expect_true(d$stationary)
  expect_identical(d$stop_reason, "all_parameters_fixed")
  expect_true(is.na(d$pi_fw_gap))
  expect_true(is.na(d$eta_kkt_residual))
  expect_false(any(d$boundary$estimated))
  expect_identical(d$information$status, "no_free_parameters")
})

test_that("analytic scores agree with finite differences including MAP penalty", {
  fit <- convergence_fit(max_iter = 2, dirichlet_alpha = c(3, 2, 2, 2))
  d <- mirage_convergence(fit); g <- d$scores; p <- fit$parameters
  objective <- function(p) {
    z <- compute_joint_quantities(
      compute_log_gene_bf(fit$variant_info$trait1, fit$posterior$Gene, p$eta_trait1),
      compute_log_gene_bf(fit$variant_info$trait2, fit$posterior$Gene, p$eta_trait2), p$pi)
    z$loglik + sum((c(3, 2, 2, 2) - 1) * log(p$pi))
  }
  h <- 1e-6
  for (j in 1:3) {
    pp <- pm <- p; pp$pi[j] <- pp$pi[j] + h; pp$pi[4] <- pp$pi[4] - h
    pm$pi[j] <- pm$pi[j] - h; pm$pi[4] <- pm$pi[4] + h
    expect_equal(unname(g$pi[j] - g$pi[4]), (objective(pp) - objective(pm)) / (2*h), tolerance = 1e-6)
  }
  for (key in c("eta_trait1", "eta_trait2")) {
    pp <- pm <- p; pp[[key]] <- pp[[key]] + h; pm[[key]] <- pm[[key]] - h
    expect_equal(unname(g[[key]]), (objective(pp) - objective(pm)) / (2*h), tolerance = 1e-6)
  }
  expect_identical(d$objective, "log_posterior_up_to_constant")
  expect_true(d$objective_non_decreasing)
  expect_true(all(diff(fit$traces$objective) >= -1e-8))
})

test_that("uninformative data have zero curvature rather than false precision", {
  d <- convergence_data(rep(0, 4))
  fit <- mirage_two_trait(d, d, 100, 100, 100, 100,
    log_bf_col_trait1 = "log_bf", log_bf_col_trait2 = "log_bf")
  z <- mirage_convergence(fit, hessian = TRUE)
  expect_true(z$converged)
  expect_true(z$stationary)
  expect_equal(z$pi_fw_gap, 0, tolerance = 1e-10)
  expect_identical(z$information$status, "computed")
  expect_false(z$information$positive_definite)
  expect_equal(z$information$numerical_rank, 0)
  expect_true(is.infinite(z$information$condition_number))
})

test_that("information matrix agrees with independent analytic pi curvature", {
  fit <- convergence_fit(estimate_pi = TRUE, fixed_eta_trait1 = .3,
    fixed_eta_trait2 = .25, max_iter = 1)
  z <- mirage_convergence(fit, hessian = TRUE)
  expect_identical(z$information$status, "computed")
  p <- fit$parameters$pi; pivot <- which.max(p); free <- setdiff(1:4, pivot)
  lb1 <- fit$gene_bf$log_B_trait1; lb2 <- fit$gene_bf$log_B_trait2
  F <- exp(cbind(0, lb1, lb2, lb1 + lb2)); den <- drop(F %*% p)
  X <- (F[, free, drop = FALSE] - F[, pivot]) / den
  expect_equal(unname(z$information$matrix), unname(crossprod(X)), tolerance = 1e-4)
})

test_that("diagnostic calls do not mutate fits and invalid controls are rejected", {
  fit <- convergence_fit(max_iter = 2); before <- serialize(fit, NULL)
  invisible(mirage_convergence(fit, hessian = TRUE))
  expect_identical(serialize(fit, NULL), before)
  expect_error(convergence_fit(max_iter = NA_real_), "positive integer")
  expect_error(convergence_fit(tolerance_loglik = 0), "finite and positive")
  expect_error(mirage_convergence(fit, stationarity_tolerance = 0), "Invalid diagnostic")
  fit$diagnostics$settings <- NULL
  expect_error(mirage_convergence(fit), "older fit")
})

test_that("diagnostics leave the previous EM trajectory and posterior unchanged", {
  old <- readRDS(testthat::test_path("fixtures", "convergence_reference.rds"))
  fit <- convergence_fit(max_iter = 25)
  expect_identical(fit$parameters, old$parameters)
  expect_identical(fit$posterior, old$posterior)
  for (key in c("loglik", "pi", "eta_trait1", "eta_trait2"))
    expect_identical(fit$traces[[key]], old[[key]])
})

test_that("fixed zero pi does not invalidate interior eta curvature", {
  fit <- convergence_fit(fixed_pi = c(.5, .5, 0, 0), max_iter = 1)
  z <- mirage_convergence(fit, hessian = TRUE)
  expect_identical(z$information$status, "computed")
  expect_true(all(is.finite(z$information$matrix)))
})
