# Current manuscript result summaries

This snapshot contains existing aggregate results only. No experiments were
rerun to create it. It supersedes older pilot summaries; historical outputs
remain unchanged. The new convergence checks are software diagnostics, not a
retroactive claim that these archived fits passed new KKT or Hessian checks.

## Simulation

500 genes, 100 variants per gene per trait; 5% baseline-frequency multiplier
indicators and fivefold multiplication; N=1000, 3000, 5000 cases **and** controls
per trait; 100 paired replicates. MIRAGE and mMIRAGE use Gamma-integrated
variant BFs and five random EM starts, selected by final log likelihood.

`simulation/method_summary.csv`, `paired_auc_differences.csv`,
`parameter_recovery.csv`, and `convergence_selected.csv` are copied verbatim
from the current manuscript's summary tables. Selected joint fits met their
EM stopping criteria in 98/100, 100/100, and 100/100 replicates, respectively.
Finite nonconverged fits were retained in performance summaries.

mTADA original rates are pre-multiplier q; oracle rates are count-generating
q_count. Both mTADA versions received true marginal risk proportions (0.05,
0.02). At N=5000 original-rate mTADA has 99 available fits. There were 189
variational diagnostic flags among 600 mTADA fits. These are limitations of the
comparison, not grounds for selecting a favorable configuration.

Mean joint empirical FDP at nominal 5% Bayesian FDR remained above nominal
for several targets: N=3000 gives 9.50%, 9.08%, 10.82% for trait 1, trait 2
and shared risk. High AUC and numerical convergence are not calibration.

Source: local `PAPER_CURRENT/Simulation/Tables/`; original experiment
`explorations/mirage_single_style_vs_mtada/random_start_sample_sizes_20261105/`.
Raw simulation datasets and bulk fitted objects are not duplicated here.

## Real-data summary

Existing continued/convergence-checked fits: single ASD 13, SCZ 286; joint ASD
19, SCZ 295; shared-risk posterior discoveries 8, at nominal 5% Bayesian FDR.
The intersection of the two marginal joint discovery lists is 11, a different
criterion from the shared-risk posterior list. The SCZ marginal gene-risk
proportion was fixed at 5%, not freely identified. ASD annotation probabilities
approached boundaries. Numerical convergence is not parameter identification
or posterior calibration.

In `joint_parameters.csv`, `strict_sensitivity` is the solution reported in
the current manuscript; `old300` is retained only as a historical comparison
and is not converged. Other rows describe initialization/tolerance checks.

Included are aggregate parameters, overlap counts, mapping counts, joint and single-trait
convergence and GO term-count summaries. No participant records, variant-level
study data, full gene posterior files or manuscript drafts are added here.

GO enrichment uses the tested 19,030-gene background, mapping to 18,560 unique
Entrez IDs. ASD BP/CC/MF significant-term counts: 76/3/15; SCZ: 17/18/12;
shared-risk: 18/13/0. The 295 SCZ Ensembl genes mapped to 297 unique Entrez IDs.

**Unresolved input limitation:** a local reproducibility audit found differing
ASC/SCHEMA genome builds and a coordinate join without liftover in the original
input preparation. Convergence checks do not resolve this issue. These real-data
outputs remain exploratory; the input harmonization must be repaired and the
analysis rerun before treating the result as a validated scientific release.

Source: local `PAPER_CURRENT/Real_Data/Tables/`; audits under
`explorations/real_data_convergence_20261104/`,
`explorations/single_real_convergence_20261104/`,
`explorations/converged_real_go_20261104/`, and
`explorations/real_data_reproducibility_20261105/`. These source directories
are retained locally and not all are distributed in this software update.
