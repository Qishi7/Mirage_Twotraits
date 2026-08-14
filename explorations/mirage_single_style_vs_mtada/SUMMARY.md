# Current results

This folder now contains two runs.

## Control-informative 30-replicate run

Output:

- `results_control_informative_30rep/`
- completed replicates: 30
- BF mode: `fast`
- genes per replicate: 500
- variants per gene: 15
- background variant probability: 0.18
- background count multiplier: 8

This is the MIRAGE-favorable design. Some noncausal variants have inflated
background counts in both cases and controls, while mTADA receives only the
baseline mutation/rate input after aggregation. MIRAGE can use the control counts
to down-weight those variants.

Mean AUC:

| target | mirage_joint | separate_mirage | mtada_aggregated |
|---|---:|---:|---:|
| trait1 | 0.747 | 0.741 | 0.605 |
| trait2 | 0.688 | 0.666 | 0.580 |
| pleiotropy | 0.790 | 0.778 | 0.593 |

Mean PP >= 0.8 discoveries and FDP:

| target | method | discoveries | observed FDP |
|---|---|---:|---:|
| trait1 | mirage_joint | 5.07 | 0.084 |
| trait1 | separate_mirage | 4.50 | 0.047 |
| trait1 | mtada_aggregated | 270.50 | 0.939 |
| trait2 | mirage_joint | 1.80 | 0.128 |
| trait2 | separate_mirage | 1.23 | 0.033 |
| trait2 | mtada_aggregated | 251.27 | 0.976 |
| pleiotropy | mirage_joint | 0.63 | 0.106 |
| pleiotropy | separate_mirage | 0.03 | 0.000 |
| pleiotropy | mtada_aggregated | 152.30 | 0.988 |

Selected pi learning:

| method | pi11 truth | mean pi11 estimate |
|---|---:|---:|
| mirage_joint | 0.010 | 0.023 |
| separate_mirage | 0.010 | 0.001 |
| mtada_aggregated | 0.010 | 0.009 |

Interpretation: under a design where controls are informative for distinguishing
true risk variants from high-background variants, MIRAGE joint clearly beats the
aggregated mTADA comparator in ranking. mTADA's pi11 estimate remains near the
truth because its marginal risk probabilities are fixed, but its posterior
discoveries are badly miscalibrated.

## Fast 30-replicate run

Output:

- `results_fast_100rep/`
- actual completed replicates: 30
- BF mode: `fast`
- genes per replicate: 500
- variants per gene: 15

Mean AUC:

| target | mirage_joint | separate_mirage | mtada_aggregated |
|---|---:|---:|---:|
| trait1 | 0.712 | 0.704 | 0.736 |
| trait2 | 0.671 | 0.629 | 0.720 |
| pleiotropy | 0.749 | 0.732 | 0.791 |

Mean discoveries:

| target | metric | mirage_joint | separate_mirage | mtada_aggregated |
|---|---|---:|---:|---:|
| trait1 | PP >= 0.8 | 4.67 | 3.47 | 5.23 |
| trait2 | PP >= 0.8 | 1.83 | 0.90 | 1.50 |
| pleiotropy | PP >= 0.8 | 1.17 | 0.03 | 0.27 |
| trait1 | bFDR 0.05 | 3.93 | 3.00 | 5.10 |
| trait2 | bFDR 0.05 | 1.07 | 0.57 | 1.33 |
| pleiotropy | bFDR 0.05 | 0.50 | 0.00 | 0.17 |

Selected pi learning:

| method | pi11 truth | mean pi11 estimate |
|---|---:|---:|
| mirage_joint | 0.010 | 0.048 |
| separate_mirage | 0.010 | 0.001 |
| mtada_aggregated | 0.010 | 0.011 |

## Integrated BF sanity check

Output:

- `results_integrated_sanity/`
- completed replicates: 3
- BF mode: `integrated`
- genes per replicate: 200
- variants per gene: 8

Mean AUC:

| target | mirage_joint | separate_mirage | mtada_aggregated |
|---|---:|---:|---:|
| trait1 | 0.677 | 0.734 | 0.676 |
| trait2 | 0.480 | 0.460 | 0.558 |
| pleiotropy | 0.517 | 0.451 | 0.591 |

Selected pi learning:

| method | pi11 truth | mean pi11 estimate |
|---|---:|---:|
| mirage_joint | 0.010 | 0.041 |
| separate_mirage | 0.010 | 0.001 |
| mtada_aggregated | 0.010 | 0.008 |

## Interpretation

This MIRAGE-single-style simulation is closer to MIRAGE than the original
mTADA-style count simulation, because the raw data are variant-level case/control
counts with category-specific risk-variant probabilities. However, the
`mtada_aggregated` comparator still receives an aggregated burden-count version
of the data, which remains favorable for mTADA-like ranking.

The current pattern is:

- mTADA aggregated has higher mean AUC in the fast 30-replicate run.
- MIRAGE joint is usually better than separate MIRAGE, especially for trait2 and
  pleiotropy discoveries.
- MIRAGE joint overestimates `pi11`; this is the main parameter-learning issue.
- mTADA's `pi11` estimate looks closer partly because the implementation fixes
  marginal risk probabilities and estimates only the overlap.

The next useful refinements are:

- add a regularized/MAP `pi` setting for MIRAGE joint;
- run a larger integrated-BF experiment after optimizing BF caching;
- compare MIRAGE joint primarily against separate MIRAGE, with mTADA aggregated
  presented as a non-native stress-test comparator.
