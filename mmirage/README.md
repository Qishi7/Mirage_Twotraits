# mMIRAGE 0.2.0

Two-trait conditional MIRAGE for gene discovery from variant-level case-control
allele counts. State order: `00`, `10`, `01`, `11`.

## Install

From the repository root, run `R CMD INSTALL mmirage`, or in R:

```r
remotes::install_github("Qishi7/Mirage_Twotraits", subdir = "mmirage")
```

## Fit from allele counts and inspect convergence

```r
library(mmirage)
# Illustrative input only; replace with your own variant-level data.
x <- data.frame(Gene = c("G1", "G1", "G2", "G2"),
                ID = paste0("v", 1:4), No.case = c(8, 3, 0, 1),
                No.contr = c(1, 1, 2, 1), category = c("C1", "C2", "C1", "C2"))
fit <- mirage_two_trait(x, x, 3000, 3000, 3000, 3000,
                       gamma_trait1 = 3, sigma_trait1 = 1,
                       gamma_trait2 = 3, sigma_trait2 = 1,
                       max_iter = 10000)
print(fit)
fit$diagnostics$convergence
head(fit$posterior)

# Optional: numerical local curvature; potentially expensive for large datasets.
checks <- mirage_convergence(fit, hessian = TRUE)
checks$information
checks$boundary
checks$warnings
matplot(fit$traces$pi, type = "l", lty = 1, xlab = "EM iteration", ylab = "pi")
```

Fitting from counts computes Gamma-integrated BFs. Supplying `log_bf_col_trait1`
or `log_bf_col_trait2` uses those BFs unchanged, so the caller must verify their
provenance. The fit does not automatically run multiple starts. When using
multiple starts, retain the highest final likelihood and report convergence for
the selected fit and all starts; do not select by AUC.

`converged` means both EM stopping tolerances passed. `stationary` separately
checks constrained first-order conditions. Boundary and curvature flags assess
additional numerical concerns; none establishes global optimality, parameter
identifiability, calibrated posterior probabilities, or FDR control.

The optional `mirage_refine()` optimizer preserves the original EM object and
returns a separate ML result; it is not automatically enabled.

See [the convergence guide](../documentation/convergence.md) and
[current manuscript result summaries](../results/current/README.md).
The script `inst/scripts/reproduce_one_case.R` is a **legacy illustrative
experiment** (15 variants/gene, fast BFs, truth-informed starting/prior values,
100 iterations), not a reproduction of the current 100-variant, random-start,
Gamma-integrated manuscript comparison. It is retained unchanged for provenance.
