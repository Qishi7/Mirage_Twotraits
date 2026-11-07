\name{mirage_refine}
\alias{mirage_refine}
\alias{mirage_simplex_optimize}
\title{Optional simplex-based maximum-likelihood refinement}
\description{Refines the unchanged conditional mMIRAGE likelihood using stored
variant Bayes factors. The original EM fit is retained as a candidate and is
never modified. No penalties or priors are added.}
\usage{
mirage_refine(fit, estimate_pi, estimate_eta, starts = list(),
  tolerance_pi = 1e-7, tolerance_eta = 1e-6, max_iter = 3000L)

mirage_simplex_optimize(log_components, pi, tolerance = 1e-7,
  max_iter = 2000L)
}
\arguments{
\item{fit}{An unpenalized \code{mirage_two_trait} EM fit. MAP fits are rejected.}
\item{estimate_pi, estimate_eta}{Explicit logical flags. Use the same fixed/free
parameter assumptions as the original analysis.}
\item{starts}{Additional feasible lists with \code{pi} (state order 00,10,01,11)
and \code{eta} (trait 1 categories followed by trait 2 categories, in the order
of the original fit's named eta vectors). Fixed parameters must be unchanged.}
\item{tolerance_pi, tolerance}{Absolute, unnormalized log-likelihood
Frank--Wolfe gap tolerance for the simplex subproblem.}
\item{tolerance_eta}{Absolute unnormalized eta box-KKT residual tolerance.}
\item{max_iter}{Iteration limit for the outer optimizer or simplex solver.}
\item{log_components}{Finite gene-by-four matrix of log component BFs,
usually columns 0, log B1, log B2, log B1 + log B2.}
\item{pi}{Feasible starting four-state simplex weights, including zeros.}
}
\details{The fixed-eta problem is concave in pi. The simplex solver uses
active-face Newton steps and pairwise vertex exchanges, allowing zero weights
to re-enter. Its certificate is max(g) - sum(pi*g), where
g[s] = sum_i f[i,s]/D[i]. This bounds the fixed-eta objective suboptimality.
It does not certify the nonconvex joint problem.

When eta is free, box-constrained L-BFGS-B profiles out pi. Starts include the
existing EM solution and two generic deterministic interior configurations.
The highest feasible log-likelihood among the original EM, all starts, and
all evaluated candidates is retained. Optimizer return codes alone do not
establish stationarity. Inspect \code{fw_gap}, \code{eta_kkt},
\code{stationary}, and the complete \code{starts} records. Failed attempts
are recorded; the original EM remains available. Extremely large component
dynamic ranges may exceed double precision and are reported as errors.
}
\value{\code{mirage_refine} returns a separate \code{mirage_optimization} list
containing parameters, posterior probabilities, BFDR tables, log-likelihood,
stationarity diagnostics and all start records. It is not a replacement EM
object and carries no EM convergence claim. \code{mirage_simplex_optimize}
returns pi, log-likelihood, the gradient, gap, certificate, termination reason,
iteration count, and objective trace.}
