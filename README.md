# mMIRAGE

An R package for joint rare-variant gene discovery from two traits using
variant-level case-control allele counts.

mMIRAGE combines trait-specific conditional MIRAGE likelihoods through four
gene-risk states: neither trait (`00`), trait 1 (`10`), trait 2 (`01`), and both
traits (`11`). Shared genes do not require identical causal variants in both
studies. Model parameters are estimated by EM with Gamma-integrated variant
Bayes factors.

## Installation

```r
remotes::install_github("Qishi7/Mirage_Twotraits")
```

Or clone the repository and run `R CMD INSTALL .`.

## Usage

Each input is a data frame with `Gene`, `ID`, `No.case`, `No.contr`, and
`category` columns. Provide case and control sample sizes for each trait and
scientifically justified Gamma relative-risk parameters.

```r
library(mmirage)
fit <- mirage_two_trait(
  trait1, trait2,
  n1_trait1 = 3000, n0_trait1 = 3000,
  n1_trait2 = 3000, n0_trait2 = 3000,
  gamma_trait1 = 3, sigma_trait1 = 1,
  gamma_trait2 = 3, sigma_trait2 = 1,
  max_iter = 10000
)
print(fit)
head(fit$posterior)
fit$bfdr$trait1
fit$diagnostics$convergence
```

A [runnable example](examples/fit_case_control.R) illustrates the input format.
Supplied log-BF columns are used unchanged; omitting them computes integrated
BFs from counts. Gamma is parameterized by mean and rate (shape = mean × rate).
The example hyperparameters are illustrative, not estimates for a new study.

## Convergence diagnostics

The fit reports the stopping reason, iteration count, final parameter and
relative likelihood changes, objective monotonicity, parameter boundaries, and
constrained first-order stationarity. Optional local curvature is available:

```r
checks <- mirage_convergence(fit, hessian = TRUE)
checks$information
checks$warnings
```

Small EM updates do not establish stationarity, parameter identification or
FDR calibration. See the [convergence guide](documentation/convergence.md) for
interpretation. The optional `mirage_refine()` optimizer returns a separate ML
result and never replaces the original EM fit automatically.

## Current analyses

- [Simulation](documentation/simulation.md): 500 genes, 100 variants/gene/trait,
  N=1000/3000/5000 cases and controls, 100 paired replicates, five random EM starts.
- [Real data](documentation/real_data.md): exploratory autism and schizophrenia
  analysis, with numerical convergence and input limitations stated explicitly.
- [Results and figures](results/README.md): current aggregate tables and the
  three manuscript figures.
- [Project website](https://qishi7.github.io/Mirage_Twotraits/).

## Repository structure

```text
R/               Inference, Bayes factors, diagnostics and evaluation utilities
man/             R function reference
examples/        Runnable input and fitting example
tests/           Software tests and small synthetic fixtures
documentation/   Model, input requirements and diagnostic interpretation
results/         Current manuscript summary tables and figures
```

## Testing

```r
# From the repository root:
pkgload::load_all()
testthat::test_local()
```

Package validation: `R CMD build .`, followed by
`R CMD check --no-manual mmirage_0.3.0.tar.gz`.

Raw study inputs and generated analysis runs are not distributed in this
software repository. License: MIT.
