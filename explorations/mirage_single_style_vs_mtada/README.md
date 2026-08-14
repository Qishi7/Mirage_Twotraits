# MIRAGE-single-style simulation vs mTADA exploration

This folder contains a small exploratory benchmark that asks what happens if the
two-trait model is tested on a variant-level simulation closer to the original
single-trait MIRAGE simulation design.

The simulation differs from the mTADA-style benchmark:

- rows are variants, not gene-category counts;
- each gene has multiple variants;
- variants have annotation categories;
- risk genes contain risk variants according to category-specific eta values;
- case and control allele counts are generated for each variant.

Simulation designs:

- `paper_like`: the earlier MIRAGE-single-style design. Case/control counts use
  the same baseline variant rate, and risk variants increase case counts.
- `control_informative`: a MIRAGE-favorable design. Some variants have inflated
  background counts in both cases and controls, while mTADA still receives the
  baseline rate. This makes control counts informative for down-weighting
  noncausal high-background variants.

Methods compared:

- `mirage_joint`: full two-trait MIRAGE fit with pi and eta estimated;
- `separate_mirage`: independent/separate MIRAGE-style baseline with fixed
  independent four-state prior;
- `mtada_aggregated`: exploratory mTADA baseline after aggregating variant-level
  case counts and allele-frequency rates to gene-category counts.

Important caveat: `mtada_aggregated` is not a native analysis for this data
generating mechanism. It ignores controls and treats aggregated case counts as
de novo-style counts, so it is included only as a stress-test comparator.

The script supports two MIRAGE variant-BF modes:

- `fast`: a binomial approximation using the category relative-risk mean;
- `integrated`: the Gamma-integrated MIRAGE variant BF from `R/two_trait_mirage.R`.

The fast mode is intended for many-replicate exploration. The integrated mode is
slower and is intended as a smaller sanity check closer to the MIRAGE model.

Run:

```sh
Rscript --vanilla explorations/mirage_single_style_vs_mtada/run_exploration.R
```

Environment variables:

- `MIRAGE_SINGLE_STYLE_REPS`, default `10`
- `MIRAGE_SINGLE_STYLE_GENES`, default `500`
- `MIRAGE_SINGLE_STYLE_VARIANTS`, default `15`
- `MIRAGE_SINGLE_STYLE_MAX_ITER`, default `100`
- `MIRAGE_SINGLE_STYLE_BF_MODE`, default `fast`, choices `fast` or `integrated`
- `MIRAGE_SINGLE_STYLE_DESIGN`, default `paper_like`, choices `paper_like` or
  `control_informative`
- `MIRAGE_BACKGROUND_PROB`, default `0.18`
- `MIRAGE_BACKGROUND_MULTIPLIER`, default `8`
- `MIRAGE_SINGLE_STYLE_OUTPUT`, default `explorations/mirage_single_style_vs_mtada/results`
- `MIRAGE_RUN_ORIGINAL_MTADA`, default `true`
