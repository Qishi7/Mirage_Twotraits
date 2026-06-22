source("R/two_trait_mirage.R")
source("R/mtada_style_simulation.R")

# Paper-profile pilot: the same four-state DNM generation scheme, two mutation
# categories, and ASD/SCZ trio counts described by Nguyen et al. (2020).
# Mutation rates are synthetic and no real study data are used.
sim <- simulate_mtada_style(
  n_genes = 19538,
  pi_trait1 = 0.05,
  pi_trait2 = 0.03,
  pi_pleiotropic = 0.015,
  n_trios_trait1 = 5122,
  n_trios_trait2 = 1077,
  gamma_mean_trait1 = c("MiD" = 20, "LoF" = 50),
  gamma_mean_trait2 = c("MiD" = 12, "LoF" = 2),
  seed = 20260619
)

result <- fit_mtada_style_mirage(sim)

print(result$metrics$auc)
print(result$metrics$pp_threshold)
print(result$metrics$bfdr)
print(result$metrics$pi)

output_dir <- Sys.getenv("MIRAGE_BENCHMARK_OUTPUT", unset = "output")
if (nzchar(output_dir) && tolower(output_dir) != "none") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(
    result$metrics$auc,
    file.path(output_dir, "mtada_style_pilot_auc.csv"),
    row.names = FALSE
  )
  write.csv(
    result$metrics$pp_threshold,
    file.path(output_dir, "mtada_style_pilot_pp08.csv"),
    row.names = FALSE
  )
  write.csv(
    result$metrics$pi,
    file.path(output_dir, "mtada_style_pilot_pi.csv"),
    row.names = FALSE
  )
  write.csv(
    result$metrics$bfdr,
    file.path(output_dir, "mtada_style_pilot_bfdr.csv"),
    row.names = FALSE
  )
}
