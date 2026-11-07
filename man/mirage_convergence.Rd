\name{mirage_convergence}
\alias{mirage_convergence}
\title{Inspect EM termination, constrained stationarity and local curvature}
\usage{mirage_convergence(fit, hessian = FALSE, stationarity_tolerance = 1e-6)}
\arguments{
\item{fit}{A \code{mirage_two_trait} fit produced by version 0.2.0 or later.}
\item{hessian}{Compute a numerical negative objective Hessian for free parameters at an interior solution. Default FALSE to avoid expensive extra likelihood evaluations.}
\item{stationarity_tolerance}{Positive tolerance for first-order residuals divided by the number of genes. This diagnostic tolerance does not alter EM fitting.}
}
\description{EM fits automatically include the basic checks in \code{fit$diagnostics$convergence}. This accessor recalculates the checks without changing parameters, posterior probabilities, or variant Bayes factors.}
\details{
\code{converged} retains the original EM definition: both maximum absolute parameter change and relative observed log-likelihood change are below their requested tolerances. The relative change is abs(new-old)/(1+abs(old)). A fit reaching its iteration limit is not automatically converged. With all parameters fixed, the reason is \code{all_parameters_fixed}, which is evaluation rather than parameter estimation.

For estimated pi the simplex Frank--Wolfe gap is max(g)-sum(pi*g), with g the objective gradient. For estimated eta the box KKT residual suppresses only outward gradients at exact endpoints. Fixed parameters are excluded. Stationarity checks use residuals per gene. Small EM updates can occur at nonstationary absorbing boundaries; the convergence and stationarity flags are intentionally separate.

The monitored objective is the observed log likelihood for ML or the existing log posterior up to a constant for \code{dirichlet_alpha} MAP fits. Both likelihood and objective traces are retained. The patch does not change either objective or the original EM stopping rule.

The optional information matrix is the negative numerical Hessian of this objective in native eta probabilities and three independent simplex coordinates (the largest pi is eliminated). For MAP fits it is posterior-objective curvature, not likelihood-only information. Boundary fits are reported as nonregular without computing an unconstrained Hessian. Eigenvalues, numerical rank and conditioning depend on parameterization, finite-difference step and numerical precision. No automatic confidence intervals are returned. Neither stationarity nor positive curvature establishes global optimality, identifiability, posterior calibration or FDR control. A regression hat matrix is not defined for this mixture-model diagnostic.
}
\value{A list with stopping status and reason, iterations, final log likelihood, final parameter and relative likelihood changes, tolerances, objective monotonicity, stationarity and raw and per-gene residuals, parameter boundary table, optional information matrix and eigenvalues, and explanatory warning strings. Warnings are returned and printed with the fit, rather than emitted on every iteration. Older saved fits without estimation metadata require refitting.}
