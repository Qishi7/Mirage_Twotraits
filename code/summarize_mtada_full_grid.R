source("R/two_trait_mirage.R")
source("R/mtada_style_simulation.R")
source("R/mtada_full_grid.R")

output_dir <- Sys.getenv(
  "MIRAGE_GRID_OUTPUT", "simulation_results/mtada_full_grid"
)
result <- collect_mtada_grid_results(output_dir)

result$metrics$threshold_group <- ifelse(
  is.na(result$metrics$threshold),
  "NA",
  format(result$metrics$threshold, scientific = FALSE, trim = TRUE)
)
summary <- aggregate(
  value ~ method + target + metric + threshold_group +
    pi_pleiotropic + rr_setting + n_trios_trait1 + n_trios_trait2,
  data = result$metrics,
  FUN = function(x) c(
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    n = sum(is.finite(x))
  )
)
summary <- cbind(
  summary[setdiff(names(summary), "value")],
  as.data.frame(summary$value)
)
summary$threshold <- suppressWarnings(as.numeric(summary$threshold_group))
summary$threshold_group <- NULL

write.csv(result$metrics, file.path(output_dir, "all_metrics.csv"),
          row.names = FALSE)
write.csv(result$task_status, file.path(output_dir, "task_status.csv"),
          row.names = FALSE)
write.csv(summary, file.path(output_dir, "summary_metrics.csv"),
          row.names = FALSE)

cat("Completed task files:", length(unique(result$metrics$task_id)), "\n")
print(summary)
