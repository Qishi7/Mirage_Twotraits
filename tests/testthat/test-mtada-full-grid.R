test_that("paper grid reproduces all Figure 2 combinations", {
  grid <- mtada_paper_grid(n_replicates = 2)
  expect_equal(nrow(grid), 2 * 4 * 9 * 2)
  expect_equal(sort(unique(grid$pi_pleiotropic)), c(0, 0.02))
  expect_equal(sort(unique(grid$n_trios_trait1)), c(1000, 2000, 5000))
  expect_equal(sort(unique(grid$n_trios_trait2)), c(1000, 2000, 5000))
  expect_equal(
    sort(unique(grid$rr_setting)),
    sort(c("high_low", "high_high", "low_low", "low_high"))
  )
})

test_that("full-grid task evaluates joint and separate methods", {
  set.seed(40)
  rates <- matrix(
    rlnorm(400, log(1e-5), 0.5),
    ncol = 2,
    dimnames = list(NULL, c("MiD", "LoF"))
  )
  task <- mtada_paper_grid(n_replicates = 1)[1, , drop = FALSE]
  result <- run_mtada_grid_task(
    task,
    mutation_rates = rates,
    run_original_mtada = FALSE
  )
  expect_true(all(c("mirage_joint", "separate") %in%
                    unique(result$metrics$method)))
  expect_true(all(c(
    "auc", "pp08_power", "pp08_fdp", "bfdr_observed_fdp",
    "bfdr_estimate", "pi11_estimate", "runtime_seconds"
  ) %in% unique(result$metrics$metric)))
})

test_that("specialized joint fit matches the four-state likelihood", {
  log_bf1 <- c(-1, 0, 2, 4)
  log_bf2 <- c(-2, 1, 0, 3)
  fit <- fit_four_state_log_bf(log_bf1, log_bf2)
  expected <- compute_joint_quantities(log_bf1, log_bf2, fit$pi)
  expect_equal(fit$tau, expected$tau, tolerance = 1e-12)
  expect_equal(sum(fit$pi), 1, tolerance = 1e-12)
  expect_true(all(diff(fit$trace_loglik) >= -1e-8))
})
