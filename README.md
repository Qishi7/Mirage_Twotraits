# Mirage_Twotraits

This repository hosts the MIRAGE Two Trait Extension: a two-trait extension of MIRAGE for joint gene-level rare-variant analysis. The motivating first use case is autism as trait 1 and schizophrenia as trait 2.

Website: https://qishi7.github.io/Mirage_Twotraits/

Repository: https://github.com/Qishi7/Mirage_Twotraits

## Model Direction

The first implementation uses a gene-level four-state model:

- `00`: neither trait
- `10`: trait 1 only
- `01`: trait 2 only
- `11`: both traits

State `11` represents gene-level pleiotropy. It does not require the same variants to affect both traits. Conditional on gene state, the two trait-specific likelihoods are independent products of trait-specific MIRAGE likelihoods.

## Repository Contents

- `R/`: first R implementation of the two-trait MIRAGE EM model and toy simulation helper.
- `R/mtada_style_simulation.R`: synthetic de novo simulation and benchmarking
  patterned after the mTADA study design.
- `code/run_mtada_style_benchmark.R`: paper-profile synthetic pilot runner.
- `code/run_mtada_full_grid.R`: resumable 7,200-dataset mTADA paper grid.
- `code/summarize_mtada_full_grid.R`: aggregate AUC, power, FDP, BFDR,
  parameter recovery, and runtime across completed tasks.
- `tests/testthat/`: unit tests and synthetic simulation smoke checks.
- `analysis/`: workflowr source pages.
- `docs/`: rendered GitHub Pages website.
- `Mirage_Twotraits.Rproj`: RStudio project file.

## Data Privacy

No real ASD, schizophrenia, WES, dbGaP, SPARK, protected, private, or large data files should be committed. Tests and examples use toy or simulated data only. Local reference documents such as papers, technical specifications, TODO files, VCF exploration reports, and implementation plans are excluded from version control.

## Running Tests

From the repository root, run:

```powershell
Rscript --vanilla tests/testthat.R
```

On the development Windows machine used for the first build, the explicit local command was:

```powershell
& 'D:\software\R-4.4.2\bin\Rscript.exe' --vanilla tests/testthat.R
```

Run the synthetic mTADA-style paper-profile pilot with:

```powershell
Rscript --vanilla code/run_mtada_style_benchmark.R
```

Set `MIRAGE_BENCHMARK_OUTPUT=none` to print results without writing CSV files.

## Full mTADA Paper Grid

The full comparison uses the public mutation-rate table from the mTADA
repository without committing that table here. Clone mTADA outside the tracked
project files, or point `MTADA_REFERENCE_DIR` to an existing local clone.

The exact Figure 2 grid contains 7,200 datasets:

- `pi_11`: 0 and 0.02
- trait-specific mRR pairs: `(105,29)/(12,2)`, `(105,29)/(87,23)`,
  `(24,6)/(12,2)`, and `(24,6)/(87,23)`
- trio counts: all nine pairs formed from 1,000, 2,000, and 5,000
- 100 repeated seeds per setting

Run a task shard on Windows:

```powershell
$env:MTADA_REFERENCE_DIR='D:\path\to\mTADA'
$env:MIRAGE_TASK_START='1'
$env:MIRAGE_TASK_END='900'
$env:MIRAGE_RUN_ORIGINAL_MTADA='true'
Rscript --vanilla code/run_mtada_full_grid.R
```

Use non-overlapping task ranges for parallel workers. Each task is saved
separately, so interrupted runs resume without repeating completed datasets.
After all shards finish:

```powershell
Rscript --vanilla code/summarize_mtada_full_grid.R
```

The original mTADA Stan model is compiled once per worker. On the development
machine, a full 19,358-gene mTADA VB fit with 5,000 iterations took about
3.5 minutes, so the original method is the dominant computational cost.
