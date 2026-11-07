# Legacy illustrative script

`reproduce_one_case.R` is retained unchanged as an early development example.
It is **not** a reproduction of the current manuscript experiment.

Its settings are 500 genes, 15 variants/gene/trait, N=3000 cases and controls,
5% background indicators and fivefold frequency multiplication, precomputed
fast BFs, truth-informed initial values, 100 EM iterations, and an independent
comparator with true fixed marginal proportions. These differ from the current
100-variant, Gamma-integrated, five-random-start comparison. In particular,
100 iterations need not converge.

To run this legacy example explicitly, install the package with
`R CMD INSTALL mmirage`, then run
`Rscript mmirage/inst/scripts/reproduce_one_case.R`. It writes to
`reproducibility_output/`; use a separate working directory to preserve any
previous outputs. For current input and diagnostic examples, see the package
README and `documentation/convergence.md` in the repository.
