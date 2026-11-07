# Manuscript results

This directory contains the current aggregate tables and publication figures.

- `simulation/`: method summaries, paired AUC differences, parameter recovery
  and selected-fit convergence for N=1000, 3000 and 5000.
- `real_data/`: final parameter estimates, final-fit convergence, gene-list
  overlap counts, GO term counts and mapping summary.
- `figures/`: the two simulation figures and real-data four-panel figure,
  each in PDF and PNG format.

See the [simulation design](../documentation/simulation.md) and
[real-data analysis](../documentation/real_data.md) for methods and limitations.
`real_data/parameter_estimates.csv` contains only the strict-tolerance fit used
in the manuscript. Convergence in archived results refers to their original
EM stopping criteria, not retrospective certification by the new diagnostics.
These are saved aggregate results; the underlying analyses were not rerun for
this repository organization.
