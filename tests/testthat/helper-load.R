repo_root <- normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_root, "R", "two_trait_mirage.R"))
source(file.path(repo_root, "R", "simulate_two_trait_mirage.R"))

toy_fit <- function(...) {
  sim <- simulate_two_trait_mirage(
    n_genes = 80,
    variants_per_gene = 6,
    pi = c("00" = 0.70, "10" = 0.10, "01" = 0.10, "11" = 0.10),
    eta_trait1 = c("1" = 0.35, "2" = 0.08),
    eta_trait2 = c("1" = 0.30, "2" = 0.06),
    signal_mean_trait1 = 3,
    signal_mean_trait2 = 3,
    null_mean = -0.8,
    seed = 1
  )
  fit <- mirage_two_trait(
    sim$trait1, sim$trait2,
    n1_trait1 = 500, n0_trait1 = 500,
    n1_trait2 = 500, n0_trait2 = 500,
    log_bf_col_trait1 = "log_bf",
    log_bf_col_trait2 = "log_bf",
    eta_init_trait1 = c("1" = 0.2, "2" = 0.05),
    eta_init_trait2 = c("1" = 0.2, "2" = 0.05),
    pi_init = c("00" = 0.75, "10" = 0.08, "01" = 0.08, "11" = 0.09),
    max_iter = 100,
    tolerance_parameter = 1e-5,
    tolerance_loglik = 1e-8,
    ...
  )
  list(sim = sim, fit = fit)
}
