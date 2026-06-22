source("R/two_trait_mirage.R")
source("R/mtada_style_simulation.R")
source("R/mtada_full_grid.R")

reference_dir <- Sys.getenv("MTADA_REFERENCE_DIR", ".tmp_mtada_reference")
mutation_rate_file <- file.path(
  reference_dir, "data", "FullDataSet_DenovoMutations_for_mTADA.txt"
)
mtada_source_file <- file.path(reference_dir, "script", "mTADA.R")
output_dir <- Sys.getenv(
  "MIRAGE_GRID_OUTPUT", "simulation_results/mtada_full_grid"
)
task_start <- as.integer(Sys.getenv("MIRAGE_TASK_START", "1"))
task_end <- as.integer(Sys.getenv("MIRAGE_TASK_END", "7200"))
run_original <- tolower(
  Sys.getenv("MIRAGE_RUN_ORIGINAL_MTADA", "true")
) %in% c("true", "1", "yes")
vb_iterations <- as.integer(Sys.getenv("MIRAGE_MTADA_VB_ITER", "5000"))

run_mtada_full_grid(
  mutation_rate_file = mutation_rate_file,
  output_dir = output_dir,
  task_ids = seq.int(task_start, task_end),
  mtada_source_file = mtada_source_file,
  run_original_mtada = run_original,
  vb_iterations = vb_iterations
)
