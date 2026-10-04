# Target-gene rank correlation results

Generated 2026-10-04T11:14:25.572780+00:00. All 55,200 competition controls are included. Gene selection is exactly the 300 entries in `pert_counts.csv`, retaining its order.

## Median per-cell Spearman correlations

| context | K562 | HCT116 | HEK293T | H1 | CD4 |
| --- | --- | --- | --- | --- | --- |
| A | 0.3229 | 0.3249 | 0.3338 | 0.3913 | 0.5129 |
| B | 0.2004 | 0.2801 | 0.3065 | 0.4625 | 0.4141 |
| C | 0.1746 | 0.2835 | 0.2602 | 0.4023 | 0.4159 |

Each correlation uses the original cell counts and the reference pooled control count vector. Only gene pairs with finite values in both vectors are ranked. Average ranks handle ties, including zeros. There is no added constant, log transform, normalization or missing-value imputation.

## Valid target-gene coverage

| reference | valid_target_genes |
| --- | --- |
| K562 | 231 |
| HCT116 | 300 |
| HEK293T | 300 |
| H1 | 300 |
| CD4 | 300 |

K562 is measured on a restricted source gene axis. Its absent targets stay missing and are excluded pairwise. Different reference pairs therefore use different gene subsets. Coverage must be considered when comparing coefficients.

## Pooled A/B/C and reference profiles

| profile | context_A | context_B | context_C | K562 | HCT116 | HEK293T | H1 | CD4 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| context_A | 1.0000 | 0.6283 | 0.5250 | 0.4563 | 0.4829 | 0.4830 | 0.5465 | 0.7468 |
| context_B | 0.6283 | 1.0000 | 0.6560 | 0.2397 | 0.3967 | 0.4252 | 0.6571 | 0.5610 |
| context_C | 0.5250 | 0.6560 | 1.0000 | 0.2561 | 0.4203 | 0.3893 | 0.6035 | 0.6014 |
| K562 | 0.4563 | 0.2397 | 0.2561 | 1.0000 | 0.6988 | 0.6993 | 0.4318 | 0.4776 |
| HCT116 | 0.4829 | 0.3967 | 0.4203 | 0.6988 | 1.0000 | 0.8307 | 0.5061 | 0.5112 |
| HEK293T | 0.4830 | 0.4252 | 0.3893 | 0.6993 | 0.8307 | 1.0000 | 0.5556 | 0.4798 |
| H1 | 0.5465 | 0.6571 | 0.6035 | 0.4318 | 0.5061 | 0.5556 | 1.0000 | 0.5764 |
| CD4 | 0.7468 | 0.5610 | 0.6014 | 0.4776 | 0.5112 | 0.4798 | 0.5764 | 1.0000 |

Pooled correlations rank summed target-gene counts across controls. They are different statistics from medians of single-cell correlations. The corresponding valid-pair counts are in `pooled_profile_valid_gene_counts.csv`. Individual A/B/C cells are not paired with cells in other contexts.

## Files and reproduction

All scripts are in the root. See [TARGET_RANK_CORRELATIONS.md](TARGET_RANK_CORRELATIONS.md) for exact commands, schemas, source scope and validation. Numeric results and provenance are in `analysis/target_rank_correlations/`. The previous broad-gene Pearson analysis remains available separately.

These are unperturbed expression associations restricted to the requested target panel. They do not measure the accuracy of CRISPRi response transfer or establish anonymous cell identities.
