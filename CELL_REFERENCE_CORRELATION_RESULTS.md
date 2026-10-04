# Control-cell reference correlation results

Generated 2026-10-04T10:54:03.622595+00:00. All **55,200 competition control cells** were analyzed: 18,400 each from A/B/C. Primary comparisons use **7,447 shared gene symbols**. No predicted cells or 2026 perturbed responses enter the calculation.

## Median per-cell Pearson correlations

| context | K562 | HCT116 | HEK293T | H1 | CD4 |
| --- | --- | --- | --- | --- | --- |
| A | 0.2938 | 0.3017 | 0.3038 | 0.3811 | 0.4821 |
| B | 0.2375 | 0.3154 | 0.2999 | 0.4635 | 0.4001 |
| C | 0.2581 | 0.3601 | 0.2969 | 0.4766 | 0.4151 |

Each value summarizes 18,400 separate cell-to-reference correlations of log1p(CPM), normalized on the same common axis. Full distributions and 5th/95th percentiles are in `context_reference_summary.csv`.

## Cells closest to each reference (%)

| context | K562 | HCT116 | HEK293T | H1 | CD4 |
| --- | --- | --- | --- | --- | --- |
| A | 0.0000 | 0.0000 | 0.0000 | 0.0054 | 99.9946 |
| B | 0.0000 | 0.0000 | 0.0000 | 99.9620 | 0.0380 |
| C | 0.0000 | 0.0000 | 0.0000 | 99.9837 | 0.0163 |

These percentages describe the largest correlation among the five candidate profiles, not known cell identities or probabilities. The per-cell CSV retains the margin between first and second place and Spearman-based choices.

## Pooled context profiles

| context | K562 | HCT116 | HEK293T | H1 | CD4 |
| --- | --- | --- | --- | --- | --- |
| A | 0.3421 | 0.4245 | 0.4288 | 0.5280 | 0.6985 |
| B | 0.2721 | 0.4613 | 0.4317 | 0.6475 | 0.5721 |
| C | 0.3284 | 0.5648 | 0.4634 | 0.7069 | 0.6313 |

Pooling counts reduces single-cell sampling zeros; pooled correlations are not averages of the per-cell correlations. A/B/C pooled profiles are also compared directly in `context_profile_correlations.csv`.

## Sensitivity checks

| context | reference | median_spearman | median_pairwise_axis_pearson | median_positive_reference_genes_pearson |
| --- | --- | --- | --- | --- |
| A | K562 | 0.3467 | 0.2915 | 0.2946 |
| A | HCT116 | 0.3387 | 0.5793 | 0.2856 |
| A | HEK293T | 0.3299 | 0.5677 | 0.2914 |
| A | H1 | 0.4166 | 0.5840 | 0.3706 |
| A | CD4 | 0.5326 | 0.6677 | 0.4693 |
| B | K562 | 0.2803 | 0.2332 | 0.2379 |
| B | HCT116 | 0.3408 | 0.5738 | 0.2973 |
| B | HEK293T | 0.3164 | 0.5509 | 0.2804 |
| B | H1 | 0.5066 | 0.6010 | 0.4573 |
| B | CD4 | 0.4613 | 0.5634 | 0.3956 |
| C | K562 | 0.2914 | 0.2551 | 0.2583 |
| C | HCT116 | 0.3782 | 0.5955 | 0.3416 |
| C | HEK293T | 0.3066 | 0.5336 | 0.2777 |
| C | H1 | 0.5155 | 0.5991 | 0.4696 |
| C | CD4 | 0.4779 | 0.5732 | 0.4060 |

Positive-reference sensitivity uses 7,309 genes with positive pooled counts in every reference. Pair-specific gene axes differ between references and cannot provide a fair closest-reference ranking.

### CD4 conditions, median per-cell Pearson

| context | reference | cells | common_genes | median_pearson | p05_pearson | p95_pearson |
| --- | --- | --- | --- | --- | --- | --- |
| A | CD4_Rest | 18400 | 7447 | 0.4701 | 0.3558 | 0.5451 |
| A | CD4_Stim8hr | 18400 | 7447 | 0.4681 | 0.3506 | 0.5439 |
| A | CD4_Stim48hr | 18400 | 7447 | 0.4806 | 0.3674 | 0.5530 |
| B | CD4_Rest | 18400 | 7447 | 0.3772 | 0.2806 | 0.4216 |
| B | CD4_Stim8hr | 18400 | 7447 | 0.3925 | 0.2892 | 0.4351 |
| B | CD4_Stim48hr | 18400 | 7447 | 0.4011 | 0.3022 | 0.4472 |
| C | CD4_Rest | 18400 | 7447 | 0.3893 | 0.3057 | 0.4588 |
| C | CD4_Stim8hr | 18400 | 7447 | 0.3933 | 0.3030 | 0.4645 |
| C | CD4_Stim48hr | 18400 | 7447 | 0.4314 | 0.3401 | 0.5081 |

## Reference scope

| reference | profile_genes | control_cells | scope |
| --- | --- | --- | --- |
| K562 | 8246 | 75094 | K562 core-control pseudobulk groups |
| HCT116 | 38584 | 891 | All passing controls in Batch1 only, not the full atlas |
| HEK293T | 38584 | 1079 | All passing controls in Batch1 only, not the full atlas |
| H1 | 18080 | 1536 | H1 2025 train release only |
| CD4 | 18129 | 34620 | CD4 non-targeting pseudobulks spanning donors, runs and all three culture conditions |

H1 uses 1,536 source controls across 48 batches. CD4 uses 384 source pseudobulks across 12 donor/condition/run strata, representing 34,620 cells. HCT116/HEK293T use first-batch controls. Thus this is a reproducible sampled-reference analysis, not a complete-atlas comparison.

## Reproduction and verification

See [CELL_REFERENCE_CORRELATIONS.md](CELL_REFERENCE_CORRELATIONS.md) for source definitions, extraction, normalization, exact commands and output schemas. `analysis/cell_reference_correlations/analysis_report.json` records source identities, parameters and hashes. Run `validate_cell_reference_correlations.py` to check complete row/identifier coverage, summaries and 45 actual-cell/reference pairs against independent NumPy/SciPy calculations; its result is saved as `validation.json`.

Basal similarity can be driven by housekeeping expression, shared zeros, sequencing depth and platform effects. It does not prove that a reference will predict CRISPRi responses accurately. Response-transfer weights require a separate held-out cell-line experiment.
