# Pairwise-valid rank correlations on the 300 target genes

## Purpose

Repeat the cell-to-reference comparison using exactly the 300 target genes in `pert_counts.csv`, with **Spearman rank correlation on original values**. Every one of the 18,400 control cells in each context A/B/C is compared with each of the five prepared external non-targeting profiles. No competition cells are subsampled.

The reference remains a pooled control-expression profile for each cell type, as in the preceding analysis. This is not an all-against-all comparison between 55,200 competition cells and individual external cells. Individual cells in A/B/C have no cross-context pairing. A separate 8-by-8 matrix compares pooled A/B/C expression with the five reference profiles.

## Method

1. Read the 300 unique symbols from `pert_counts.csv`, keeping its original order. Confirm that all are in the competition's measured gene axis.
2. Read each previously prepared reference's `control_count_sum` values, verifying the profile SHA-256 against its provenance record. Gene matching uses exact symbols.
3. Reindex each reference to the 300 target symbols. A target absent from that reference is represented as missing, not zero. K562 provides 231 targets; HCT116, HEK293T, H1 and CD4 provide all 300.
4. For each cell/reference pair, keep only genes for which **both values are finite**. Remove missing/nonfinite pairs before computing ranks.
5. Rank the remaining original values with average ranks for ties. **Measured zeros are valid and remain in the calculation.** One HCT116 target has a measured pooled zero; it is not omitted. The original competition counts are finite nonnegative integers; sparse unstored entries are measured zeros.
6. Compute the Pearson correlation of these two rank vectors (Spearman rho). This is mathematically the rank statistic, not Pearson correlation of raw expression.
7. Record the valid-pair count and status. Fewer than two valid pairs, or a constant rank vector, produce an undefined correlation (blank CSV value/NaN), with an explicit status. Undefined results are never replaced with zero. Summary quantiles exclude undefined correlations and report their counts.

**There is no addition of 1 or any other pseudocount, log transform, CPM normalization, positive-only filtering, or missing-value imputation.** K562's 69 absent target symbols are excluded pairwise. Because the gene subset differs between reference pairs, interpret coefficient comparisons with the coverage counts alongside them.

Rank correlation is invariant to multiplying either entire expression vector by a positive constant. Different library sizes therefore do not require a scaling step for this statistic. Count sparsity and many tied zeros still affect the biological interpretation.

## Reference scope

This run reuses the verified references from [CELL_REFERENCE_CORRELATIONS.md](CELL_REFERENCE_CORRELATIONS.md):

- K562: original pseudobulk raw means multiplied by control-group cell membership and pooled; 75,094 annotated control cells across 503 groups.
- HCT116: 891 passing controls from the first cached Orion batch.
- HEK293T: 1,079 passing controls from the first cached Orion batch.
- H1: 1,536 deterministic control samples, 32 per batch across 48 batches of the public 2025 train release.
- CD4: 384 non-targeting pseudobulks, spanning 12 donor/condition/run strata and representing 34,620 cells; all three culture conditions are pooled.

These are sampled or batch-limited external references. No perturbed 2026 challenge measurements are used. The analysis does not establish anonymous cell identities or predict CRISPRi response-transfer accuracy.

## Root code and outputs

- `analyze_target_rank_correlations.py`: load, align, compute, summarize, and regenerate the root results Markdown.
- `test_target_rank_correlations.py`: independent SciPy/Pandas oracles for missing values, ties, measured zeros, re-ranking after pair omission, insufficient pairs, constant vectors and positive scaling.
- `validate_target_rank_correlations.py`: check complete output coverage and compare selected actual-cell correlations and every pooled-profile matrix entry against independent SciPy calculations.

All numeric outputs are in `analysis/target_rank_correlations/`:

| File | Contents |
|---|---|
| `per_cell_target_rank_correlations.csv` | 55,200 rows, retaining context, cell ID, guide ID and original count metadata; five Spearman values, five valid-gene-pair counts and five status columns |
| `context_reference_summary.csv` | Median, mean, quantiles, defined/undefined correlation counts and valid-gene counts for each context/reference |
| `guide_reference_summary.csv` | Per-guide medians and defined-correlation counts |
| `target_reference_raw_profiles.csv` | Actual values used for each of the 300 targets, with blanks for missing reference genes |
| `target_gene_coverage.csv` | Valid-value flags for each target/reference |
| `pooled_raw_target_profiles.csv` | Summed A/B/C counts and reference profiles on the 300-target axis |
| `pooled_profile_rank_correlations.csv` | Pairwise Spearman matrix for pooled A/B/C and the five reference profiles |
| `pooled_profile_valid_gene_counts.csv` | The number of finite gene pairs for each matrix entry |
| `analysis_report.json` | Method, source scope, versions, input/output hashes, target coverage and undefined-correlation counts |
| `validation.json` | Independent complete-axis and numerical validation results |

Results tables are regenerated in [TARGET_RANK_CORRELATION_RESULTS.md](TARGET_RANK_CORRELATION_RESULTS.md). The earlier broad-gene analysis remains in its original directory.

## Reproduce in Windows Command Prompt

Use the same original competition files in the project root and the pinned environment from the previous analysis. If the external profile cache is missing, reproduce its preparation using `download_reference.py`, `pilot_orion.py` and `prepare_reference_controls.py`, as documented in [CELL_REFERENCE_CORRELATIONS.md](CELL_REFERENCE_CORRELATIONS.md).

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m unittest test_target_rank_correlations -v
.\.competition_env\Scripts\python.exe analyze_target_rank_correlations.py
.\.competition_env\Scripts\python.exe validate_target_rank_correlations.py
```

Optional computational settings:

```cmd
.\.competition_env\Scripts\python.exe analyze_target_rank_correlations.py --batch-size 1024 --min-valid-pairs 2
```

Batch size changes computation only. Raising the minimum-pair threshold changes when a correlation is defined and is recorded in provenance. The calculation generates no new prediction or leaderboard upload.

## Completed run and validation

The completed run contains 55,200 cell rows and no undefined cell/reference correlations. Median rank correlations favor CD4 for A (0.5129) and H1 for B (0.4625). C has a small per-cell median lead for CD4 (0.4159 versus H1 0.4023), while its pooled profile is nearly tied and slightly favors H1 (0.6035 versus CD4 0.6014). This difference between summaries does not establish a cell identity.

All output cell identifiers and guides, summary medians, valid-pair counts and CSV checksums passed validation. Independent SciPy calculations matched 45 cell/reference pairs read from the original H5AD counts and all 64 pooled matrix entries. Three focused unit tests passed. Detailed evidence is saved in `analysis/target_rank_correlations/validation.json` and `test_results.txt`.
