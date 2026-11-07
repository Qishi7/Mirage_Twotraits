# Model and input

For each gene the four-state mixture represents neither trait, trait 1 only,
trait 2 only or shared risk, in the order `00`, `10`, `01`, `11`. Within a risk
gene, each variant is causal with annotation-specific probability eta; within
a nonrisk gene all variants are noncausal. The two studies are conditionally
independent given gene state.

For total allele count T and relative risk gamma, the conditional case count is
Binomial(T, N_case*gamma/(N_case*gamma + N_control)). Under the null gamma=1.
Under the alternative gamma is integrated over its annotation-specific Gamma
prior, parameterized by mean and rate. Gene BFs are products of
`1 - eta + eta * variant_BF`. The gene-state marginal contribution is
`pi00 + pi10*B1 + pi01*B2 + pi11*B1*B2`.

EM estimates pi and eta. Specify fixed parameters explicitly when appropriate.
Sample sizes, Gamma hyperparameters, and supplied variant BFs are inputs rather
than parameters estimated by this EM. Multiple starts must be run explicitly;
compare final likelihoods rather than downstream ranking metrics.

## Input columns

| Column | Meaning |
|---|---|
| `Gene` | Gene identifier, consistent between the two studies |
| `ID` | Variant identifier within the study |
| `No.case` | Nonnegative integer alternate allele count in cases |
| `No.contr` | Nonnegative integer alternate allele count in controls |
| `category` | Functional annotation category |

Different column names can be supplied through the function arguments. A union
of genes is used by default. Missing study-specific variants are not filled
with artificial zero counts. A gene without variants in one study has a neutral
BF of one in that study.

Only the caller can establish that sample exposures, annotations and genome
builds are harmonized. The software does not perform liftover or study-level
variant QC. When allele numbers vary substantially, justify the fixed exposures
used by the conditional likelihood.

## Output

`posterior` contains four-state and marginal probabilities. `bfdr` contains
sorted gene probabilities and cumulative posterior false-discovery estimates.
`parameters` contains pi and eta. `traces` and `diagnostics` support numerical
assessment. See [convergence](convergence.md) before interpreting a fitted model.
