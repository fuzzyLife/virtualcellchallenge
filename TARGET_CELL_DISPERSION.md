# Cell-level correlation dispersion on the target-gene panel

## Purpose and scope

Measure variation among the 18,400 unperturbed cells in each of A/B/C using exactly the 300 targets in `pert_counts.csv`. This extends the per-cell correlations in [TARGET_RANK_CORRELATIONS.md](TARGET_RANK_CORRELATIONS.md). The external K562, HCT116, HEK293T, H1 and CD4 data remain their previously prepared pooled control profiles. Actual individual-cell comparisons in this analysis are between the supplied A/B/C cells, within and across contexts. There is no finer cell-state annotation in the supplied metadata.

## Step-by-step calculation

1. Verify the preceding correlation CSV and reference-vector hashes, and each original `context_*.h5ad` checksum. Match all cell IDs and the requested target-gene order before reading the raw counts.
2. Retain the preceding five correlations for every cell. Summarize their across-cell sample variance and SD (ddof=1), median, 5th/95th percentiles, IQR and median absolute deviation. The percentile ranges describe variation among observed cells, not confidence intervals for a mean.
3. Recompute each cell's correlations to all five external profiles on the SAME 231 target genes measured by every reference. Use these coefficients for highest-reference fractions and the highest-minus-second margin. Keep the original available-gene correlations as well (231 K562, 300 others). Missing values are omitted before ranking; measured zeros remain tied observations. No addition of 1, transformation, normalization or imputation is applied.
4. Compute each cell's correlation with its own context's raw-count pooled profile after removing that cell's counts. This leave-one-out comparison avoids counting the cell in its comparator. It is a comparison to a pooled profile, which has different noise from another individual cell.
5. Rank all 300 counts in each A/B/C cell with average ranks for ties, center the ranks, and divide by their Euclidean norm. For finite, nonconstant rows, the dot product of these unit rank vectors is precisely Spearman rho. An undefined constant row is excluded from valid peer sets. Nonfinite original count vectors are rejected: the exact moment formula requires a shared finite axis. Original counts in this run are all finite and nonconstant.
6. Compute each cell's exact mean and population SD of correlations to EVERY valid cell in each context, excluding self. For a target group with unit rank matrix Z, sum of correlations to row z is `z @ sum(Z)`, and sum of squared correlations is `z @ (Z.T @ Z) @ z`. Subtract the self-correlation and its square when the query belongs to the group. Divide by the remaining peer count; variance is second moment minus squared mean. This covers all pairs without storing a billions-entry matrix. Self-correlations never contribute. One valid peer has zero SD; small negative floating-point variance residuals are clipped to zero.
7. Aggregate the exact moments into the mean and population SD of each A-A, A-B, A-C, B-B, B-C and C-C cell-pair distribution. Same-context unordered pairs are counted once.
8. For pair-distribution percentiles and individual-pair CSVs, sample 20,000 pairs per group uniformly WITH replacement using seed 20261004. Reject self-pairs, store same-context pairs in index order, and record both IDs. Duplicates are possible and their number is reported. Sample percentiles are approximate; the means and SDs above use all valid pairs exactly.
9. For each target gene and context, calculate raw-count mean, sample variance, SD, zero fraction and Fano factor (variance/mean; undefined for zero mean). Raw-count dispersion depends on library depth and does not isolate biological variation.
10. Check Spearman associations between correlation metrics and full-library counts / number of nonzero target genes. Also summarize by non-targeting guide. These are descriptive checks, not proof of biological subtypes or batch correction.

## Outputs

Root code: `analyze_target_cell_dispersion.py`, `test_target_cell_dispersion.py`, `validate_target_cell_dispersion.py`. Root findings: [TARGET_CELL_DISPERSION_RESULTS.md](TARGET_CELL_DISPERSION_RESULTS.md).

All numeric and figure outputs are in `analysis/target_cell_dispersion/`:

| File | Contents |
|---|---|
| `per_cell_dispersion.csv` | One row for each of 55,200 cells; previous correlations, five matched-axis correlations, highest reference and margin, leave-one-out pooled correlation, and mean/SD/valid peer count for each A/B/C peer group |
| `correlation_dispersion_summary.csv` | Across-cell variance, SD, quantiles, IQR and MAD of each comparison |
| `cell_pair_dispersion_summary.csv` | Exact all-pair means/SDs/counts and sampled percentiles |
| `sampled_individual_cell_pair_correlations.csv` | 120,000 individual cell-pair correlations, IDs, contexts and valid gene counts |
| `highest_reference_fractions.csv` | Highest-reference counts and fractions on the same 231 genes, including tie and undefined categories |
| `target_gene_count_dispersion.csv` | Raw-count variance and sparsity for 900 context/gene combinations |
| `count_depth_associations.csv` | Correlation metrics versus library counts / detected target genes |
| `guide_dispersion.csv` | Per-guide median similarity, reference contrast and depth |
| `gene_axis.csv` | The 300 targets and common-reference axis flags |
| `cell_dispersion.png`, `cell_dispersion.pdf` | Four-panel scientific figure |
| `analysis_report.json` | Method, input/output hashes, seed, pair sampling, scope and undefined-vector count |
| `validation.json` | Independent validation evidence |

SD in the per-cell peer columns describes the spread of that cell's correlations to its peers. SD in the summary table describes variation of a chosen metric across cells. These have different meanings. Reference profiles are batch-limited or sampled as documented in the preceding study; similarities cannot establish the anonymous contexts' cell identities or CRISPRi transfer quality.

## Reproduce in Windows Command Prompt

Use the existing pinned environment, original competition files and verified target-profile outputs from the preceding study. If those outputs are missing, first reproduce [TARGET_RANK_CORRELATIONS.md](TARGET_RANK_CORRELATIONS.md).

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m unittest test_target_cell_dispersion -v
.\.competition_env\Scripts\python.exe analyze_target_cell_dispersion.py
.\.competition_env\Scripts\python.exe validate_target_cell_dispersion.py
```

Optional sampling settings affect sampled pairs and percentiles only, not the exact means/SDs:

```cmd
.\.competition_env\Scripts\python.exe analyze_target_cell_dispersion.py --sample-pairs 20000 --seed 20261004
```

The validator rereads original H5AD counts, checks every cell ID and prior correlation, all sampled pairs against freshly computed ranks, 36 individual pairs against SciPy `spearmanr`, 45 matched-axis reference pairs and nine leave-one-out pairs against SciPy. For nine selected cells, it explicitly constructs all peer correlations to each context and checks the saved mean and SD (27 checks). All 900 raw-count variance rows, summary statistics, pair counts and output checksums are checked.

## Completed-run interpretation

Within-context cell-pair mean Â± population SD is A 0.4698 Â± 0.0863, B 0.4949 Â± 0.0969, and C 0.4761 Â± 0.0828. Cross-context means are A-B 0.3168, A-C 0.2378 and B-C 0.3160. B has the greatest spread in individual cell-pair similarity, including a lower tail visible in the figure.

On the common 231-gene reference axis, CD4 is highest for 97.94% of A cells; H1 for 62.52% of B (CD4 36.98%); and CD4 for 68.09% of C (H1 31.88%). These fractions differ from comparisons on all 300 genes because the gene panel changes. They are not subtype proportions.

The Spearman association between total library counts and cell-to-own-context correlation is 0.9085 (A), 0.8357 (B), and 0.8881 (C). Correlations to H1/CD4 also depend on count depth (0.5604â€“0.8105). Positive scaling alone cannot change a rank correlation, but depth-associated detection and tied zeros can. These descriptive associations do not measure how much variance is technical or biological.

Validation passed for all 55,200 cell rows, 120,000 sampled pair rows and 900 gene-variance rows. Independent checks included 27 explicitly enumerated all-peer mean/SD comparisons, 45 shared-axis reference correlations, nine leave-one-out correlations and 36 individual pairs against SciPy. Five focused tests passed; evidence is in `analysis/target_cell_dispersion/validation.json` and `test_results.txt`.
