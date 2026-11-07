# mMIRAGE convergence diagnostics

Version 0.2.0 adds diagnostics to the existing EM fit without changing the
likelihood, priors, BF integration, EM updates, starting values, stopping rule
or posterior probabilities.

```r
print(fit)
d <- fit$diagnostics$convergence
# Equivalent accessor; hessian=TRUE additionally computes local curvature.
d <- mirage_convergence(fit, hessian = TRUE)
```

## How to read the output

| Output | Interpretation |
|---|---|
| `converged`, `stop_reason` | Both EM tolerances passed, iteration limit reached, or all parameters were fixed. |
| `final_parameter_change` | Maximum absolute change in pi and eta in the last update. |
| `final_relative_loglik_change` | `abs(new-old)/(1+abs(old))` for observed log likelihood, including the first iteration's change from initialization. |
| `objective_non_decreasing` | No decrease exceeding 1e-8 in the fitted objective, including initialization-to-first-update. |
| `stationary` | Per-gene simplex and box KKT residuals pass a separate tolerance (default 1e-6). |
| `pi_fw_gap`, `eta_kkt_residual` | Raw constrained first-order residuals for estimated parameters only; fixed blocks are NA. |
| `boundary` | pi and eta estimates, whether they were estimated, and proximity within 1e-8 to 0 or 1. |
| `information` | Optional negative objective Hessian, eigenvalues, numerical rank and condition number for interior free parameters. |
| `warnings` | Explanatory strings for failed tolerances, decreases, nonstationarity, boundaries or unavailable diagnostics. Also shown by `print(fit)`. |

The original top-level fields `fit$diagnostics$converged`, `n_iter`,
`loglik_non_decreasing`, and all existing parameter traces remain available.
New traces include signed and relative log-likelihood changes, the fitted
objective, and signed objective changes; starting parameters are saved under
`fit$traces$initial_parameters`. For the already-supported Dirichlet MAP option,
objective monotonicity and stationarity include the log prior, while the
original stopping rule still uses relative observed log-likelihood change.
This distinction prevents calling a legitimate MAP likelihood decrease an
objective failure.

## First-order checks

Let `F[i,s] = (1, B1, B2, B1*B2)[s]` and
`D[i] = sum_s pi[s]*F[i,s]`. For ML,
`g[s] = sum_i F[i,s]/D[i]`. The simplex Frank--Wolfe gap is
`max(g) - sum(pi*g)`, truncated at zero for rounding. For the existing proper
Dirichlet MAP option the derivative `(alpha[s]-1)/pi[s]` is added where alpha
exceeds one. The gap bounds fixed-eta objective suboptimality by concavity;
it does not certify the nonconvex joint problem.

For each eta category the likelihood derivative is the sum over its variants
of `w_gene * (BF_variant - 1)/(1 - eta + eta*BF_variant)`. It is evaluated on
the log scale, including exact endpoints. At eta=0, negative gradients point
outside the box and are suppressed; at eta=1 positive gradients are suppressed.
Interior gradients and all inward gradients remain. The maximum absolute
remaining gradient is the box KKT residual. Residuals are divided by the number
of genes for the `stationary` flag; raw values are also returned.

A zero pi or eta can be absorbing for EM even with an inward ascent direction.
Consequently `converged=TRUE` and `stationary=FALSE` is possible and meaningful.
`stop_reason=all_parameters_fixed` means evaluation, not successful estimation.

## Optional information matrix

A regression hat matrix is not a convergence diagnostic for this mixture.
Instead, `hessian=TRUE` computes negative objective curvature using three native
simplex coordinates (the largest pi is eliminated) and all free eta values.
It uses `stats::optimHess` with an interior-aware finite-difference step. Fixed
parameters are excluded. Boundary solutions return `boundary_not_regular`, and
all-fixed fits return `no_free_parameters`; no invalid boundary standard errors
are manufactured. With MAP the matrix is log-posterior curvature, explicitly
labelled, not likelihood-only observed information.

Numerical rank and positive definiteness use the threshold
`1e-7 * max(1, max(abs(eigenvalues)))`. Conditioning depends on coordinates and
finite-difference step. Treat it as a local warning rather than proof of
identification, and repeat with independent starts or further mathematical
analysis where needed. No automatic confidence intervals are supplied.

Neither small EM steps, first-order stationarity nor a positive-definite Hessian
establishes global optimality or posterior/FDR calibration. Inspect all warnings,
parameter trajectories, likelihoods and scientifically justified sensitivities.
The diagnostics do not alter or automatically refine a fit. Old serialized fits
without estimation metadata must be refitted; flags are not inferred by
executing old call expressions.

## Export a compact record

```r
d <- mirage_convergence(fit)
write.csv(data.frame(converged = d$converged, stop_reason = d$stop_reason,
                     iterations = d$n_iter, loglik = d$final_loglik,
                     parameter_change = d$final_parameter_change,
                     relative_loglik_change = d$final_relative_loglik_change,
                     stationary = d$stationary,
                     pi_fw_gap_per_gene = d$pi_fw_gap_per_gene,
                     eta_kkt_per_gene = d$eta_kkt_per_gene),
          "convergence.csv", row.names = FALSE)
write.csv(d$boundary, "parameter_boundaries.csv", row.names = FALSE)
```

Tests cover iteration exhaustion, fixed blocks, an absorbing nonstationary
boundary, finite-difference gradient checks with MAP, singular information under
uninformative BFs, analytic pi curvature, unchanged fit results, and invalid
controls. Test fixtures contain only synthetic software-test examples.
