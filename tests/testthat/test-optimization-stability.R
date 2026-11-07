stored_optimization_fit <- function() {
  x <- readRDS(testthat::test_path("fixtures", "optimization_stored.rds"))
  mirage_two_trait(x$trait1, x$trait2, 12000, 12000, 12000, 12000,
    log_bf_col_trait1 = "log_bf", log_bf_col_trait2 = "log_bf",
    max_iter = 5)
}

stored_optimization_model <- function(fit) {
  c1 <- names(fit$parameters$eta_trait1); c2 <- names(fit$parameters$eta_trait2)
  n <- nrow(fit$posterior)
  .mirage_optimization_model(
    .mirage_optimization_trait(fit$variant_info$trait1, n, c1),
    .mirage_optimization_trait(fit$variant_info$trait2, n, c2), n, length(c1), length(c2))
}

test_that("native simplex escapes a vertex and certifies the fixed-eta maximum", {
  fit <- stored_optimization_fit(); model <- stored_optimization_model(fit)
  eta <- c(fit$parameters$eta_trait1, fit$parameters$eta_trait2)
  L <- model$components(eta)$L
  # A vertex which cannot be optimal: the best pure component has greater LL.
  worst <- which.min(colSums(L)); p <- numeric(4); p[worst] <- 1
  r <- mirage_simplex_optimize(L, p)
  expect_true(r$certified)
  expect_lte(r$fw_gap, 1e-7)
  expect_gt(r$loglik, sum(L[, worst]) + 1e-5)
  expect_true(all(r$pi >= 0))
  expect_equal(sum(r$pi), 1, tolerance = 1e-14)
  # Concavity certificate independently recomputed with the original engine.
  b <- model$evaluate(r$pi, eta)
  expect_equal(r$loglik, b$ll, tolerance = 1e-11)
  expect_equal(sum(r$pi * b$grad_pi), nrow(L), tolerance = 1e-10)
  expect_lte(max(b$grad_pi) - nrow(L), 1e-7)
  for (i in 1:4) {
    v <- numeric(4); v[i] <- 1
    expect_lte(sum(L[, i]) - r$loglik, r$fw_gap + 1e-10)
  }
})

test_that("compressed likelihood and eta gradient agree with the unchanged core", {
  fit <- stored_optimization_fit(); model <- stored_optimization_model(fit)
  eta <- c(fit$parameters$eta_trait1, fit$parameters$eta_trait2)
  pi <- fit$parameters$pi
  z <- model$evaluate(pi, eta)
  expect_equal(z$ll, tail(fit$traces$loglik, 1), tolerance = 1e-10)
  for (i in seq_along(eta)) {
    ep <- em <- eta; ep[i] <- ep[i] + 1e-6; em[i] <- em[i] - 1e-6
    fd <- (model$evaluate(pi, ep)$ll - model$evaluate(pi, em)$ll) / 2e-6
    expect_equal(z$grad_eta[i], fd, tolerance = 1e-6)
  }
  for (endpoint in c(0, 1)) {
    e <- rep(endpoint, length(eta)); z <- model$evaluate(pi, e)
    c1 <- stats::setNames(e[1:3], names(fit$parameters$eta_trait1))
    c2 <- stats::setNames(e[4:6], names(fit$parameters$eta_trait2))
    original <- compute_joint_quantities(
      compute_log_gene_bf(fit$variant_info$trait1, fit$posterior$Gene, c1),
      compute_log_gene_bf(fit$variant_info$trait2, fit$posterior$Gene, c2), pi)
    expect_equal(z$ll, original$loglik, tolerance = 1e-10)
  }
})

test_that("refinement retains EM, flags constraints and never mutates fits or BFs", {
  fit <- stored_optimization_fit(); before <- serialize(fit, NULL)
  r <- mirage_refine(fit, TRUE, TRUE, max_iter = 100)
  expect_gte(r$ll, tail(fit$traces$loglik, 1) - 1e-10)
  expect_equal(serialize(fit, NULL), before)
  expect_true(all(r$pi >= 0)); expect_equal(sum(r$pi), 1, tolerance = 1e-12)
  expect_true(all(r$eta >= 0 & r$eta <= 1))
  expect_equal(rowSums(r$posterior[, 2:5]), rep(1, nrow(r$posterior)), tolerance = 1e-12)
  expect_gte(length(r$starts), 4)
  objectives <- vapply(r$starts, function(s) if (is.null(s$ll)) -Inf else s$ll, numeric(1))
  expect_equal(r$ll, max(objectives))
  fixed <- mirage_refine(fit, FALSE, FALSE)
  expect_equal(fixed$pi, unname(fit$parameters$pi), ignore_attr = TRUE)
  expect_equal(fixed$eta, unname(c(fit$parameters$eta_trait1, fit$parameters$eta_trait2)))
  expect_error(mirage_refine(fit, TRUE, FALSE,
    starts = list(list(pi = rep(.25, 4), eta = rep(.7, 6)))), "fixed eta")
  fit$call$dirichlet_alpha <- quote(c(2, 2, 2, 2))
  expect_error(mirage_refine(fit, TRUE, TRUE), "MAP")
})
