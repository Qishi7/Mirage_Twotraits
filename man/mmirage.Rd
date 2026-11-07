\name{mmirage}
\alias{mmirage-package}
\alias{add_variant_log_bf}
\alias{bfdr_table}
\alias{binary_auc}
\alias{gene_burden_test}
\alias{log_variant_bf_integrated}
\alias{mirage_two_trait}
\alias{single_trait_posterior_from_log_bf}
\title{Two-trait conditional MIRAGE model and supporting utilities}
\description{Fits a four-state gene-risk mixture from variant-level case-control allele counts or supplied log variant Bayes factors. The states are ordered 00, 10, 01, 11. Variant effects are integrated using an annotation-specific Gamma model when log Bayes factors are not supplied.}
\details{\code{mirage_two_trait} accepts two data frames with Gene, ID, No.case, No.contr and category columns, four case/control sample sizes, and optional annotation-specific gamma and sigma parameters. Fixed pi and eta arguments override estimation. The returned fit contains parameters, gene-state and marginal posterior probabilities, BFDR tables, gene Bayes factors, variant information, parameter and likelihood traces, and convergence diagnostics. Consult \code{formals(mirage_two_trait)} for all fitting options and the repository diagnostic guide for runnable examples. The two relative-risk hyperparameters specify Gamma mean and rate, with shape equal to their product.

The EM likelihood, updates and variant Bayes factors are unchanged. A supplied log-BF column is used as-is; callers must verify its provenance. The public package contains the inference software; manuscript results and methodology are documented separately in the repository.
}
\seealso{\code{\link{mirage_convergence}}, \code{\link{mirage_refine}}}
