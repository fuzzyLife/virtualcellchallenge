# Audit and clustering of the user-provided 300-target archive

## Inputs and scope

The root input `vcc_2026_300_target_matrix.zip` contains `vcc_2026_300_target_matrix.csv`. Its columns must be `context` and exactly the 300 unique `target_gene` symbols in `pert_counts.csv`. The expected rows are 18,400 unperturbed controls per context A/B/C, 55,200 cells in total. These are baseline controls, not measured CRISPRi response profiles.

The original reference counts are `context_A.h5ad`, `context_B.h5ad` and `context_C.h5ad`, each 18,400 cells by 18,533 genes. `gene_names.csv` defines their common gene axis. Original input hashes must agree with the earlier full-matrix clustering audit. The supplied ZIP remains unchanged and is treated as local input, excluded from Git; its member and archive SHA-256 checksums are tracked in the audit.

Root scripts are `analyze_target_archive.py`, `test_target_archive.py` and `validate_target_archive.py`. The analysis reuses the reviewed `rank_unit_rows` and converged spherical-clustering implementation in `cluster_full_rank_matrix.py`. Results are written under `analysis/target_archive_audit/`; the human-readable result is [TARGET_ARCHIVE_RESULTS.md](TARGET_ARCHIVE_RESULTS.md).

## 1. Validate the archive and every count

Check ZIP CRC, unambiguous CSV membership, unique column names, complete target coverage, row counts and valid context labels. All values must be finite, nonnegative and integer-valued. CSV float notation such as `2.0` represents an integer raw count, not normalized expression.

Read original matrices in 512-cell batches. Select the 300 columns by their names, and compare every value exactly with the archive in within-context row order: 55,200 x 300 = 16,560,000 comparisons. Record each original SHA-256, mismatched-value count and maximum absolute difference. Any mismatch stops analysis.

**The archive does not contain cell IDs.** Exact agreement of the complete 300-gene vector in source order verifies the ordered matrix contents. Original IDs and non-targeting guide IDs are attached from the corresponding original row; this is a verified positional mapping, not independent identification from an archived cell-ID field. `archive_row_index` is zero-based. This limitation matters if the archive is subsequently shuffled.

## 2. Cluster all cells together using only these genes

Rank each cell's 300 raw counts with average ranks for ties. Keep measured zeros. Center each rank vector and divide by its Euclidean norm; the dot product between two such vectors is Spearman correlation. No pseudocount, log transform, missing-value filling, count normalization or other gene enters clustering. Constant or nonfinite rows are rejected.

Fit MiniBatchKMeans initialization followed by spherical refinement on all 55,200 cells until assignments stop changing. Evaluate K=2,3,4,6; repeat K=3 with a second seed. Context labels are withheld from fitting. Attach A/B/C names to the primary three clusters afterward using the maximum-agreement one-to-one assignment. Report full composition, context ARI/NMI, repeat ARI and silhouettes on a fixed 1,500-cell sample using distance `1 - Spearman rho`. The sample silhouette alone does not choose a biologically true number of groups.

Fit eight PCA components to the 300-gene rank vectors for display only. Cluster fitting occurs in the original 300-dimensional rank geometry. Preserve every cell, its original ID, non-targeting guide, count depth, target detection, cluster labels, centroid similarity/margin and display coordinates. Compare discordant cells with the previous full-18,533-gene clustering.

## 3. Identify enriched genes and check depth

For marker comparison only, calculate per-cell CPM using total counts over the **full 18,533 measured genes**, not the sum of the 300 selected genes. Do not add any pseudocount or take a log before ranking. Compute mean/median raw counts and CPM, detection fractions, pairwise rank AUCs, AUC versus the other cells, and CPM mean fold relative to the higher mean in the other two groups.

Pairwise AUC includes ties at half weight. An AUC of 0.5 indicates no rank separation, and 1 indicates that the gene has greater CPM in every cell of the first group than in the other group, with no ties. The minimum of the two pairwise AUCs prevents a marker driven only by one easy comparison.

The descriptive `panel_enriched` flag requires all of:

- Minimum pairwise AUC at least 0.8 against both other groups.
- Mean CPM at least twice the larger other-group mean.
- Measured count above zero in at least 25% of the group's cells.

The separate `detection_restricted` flag requires detection in at least 50% of the group and no more than 5% in either other group. This describes these observed controls, not universal tissue specificity. Markers are ranked by minimum pairwise AUC, then mean fold and gene name. Complete 900-row tables are produced both for blind clusters and for the original contexts. A gene need not meet either threshold even when it contributes to multigene clustering.

Fold ratios with a zero denominator and a positive numerator are infinite; all-zero comparisons are undefined/NaN. Zero numerator gives a zero ratio and the log2 column is left undefined. No finite pseudocount estimate replaces these cases. The log2 column is a reporting transformation of the mean ratio, not an input to rank tests or clustering.

Split cells into five equally sized depth strata independently within each original context, using full-library depth and stable rank-based tie handling. Summarize selected top markers and the requested VIM/S100A11/FXYD5/EZH2/SMARCA5/E2F1 genes within every stratum. This checks robustness to depth without claiming to remove every technical or biological confounder.

No population-level significance claim is based on treating 18,400 cells as 18,400 independent biological replicates. Context and batch may be confounded; the selected CRISPRi target panel is not an unbiased lineage gene panel. This analysis characterizes baseline differences and does not validate perturbation prediction accuracy.

## 4. Check the biological interpretation

Use an explicit broader lineage and proliferation panel from the original 18,533-gene matrices, including CD3D/CD3E/TRAC, epithelial keratins/EPCAM/CDH1, mesenchymal-associated genes and cell-cycle genes. List any unmeasured requested markers rather than imputing them. These additional genes never enter the 300-gene clustering. Means, medians, detection fractions and full-library CPM are saved with panel membership and whether each gene belongs to the 300 targets.

The interpretation and primary-source links are in [TARGET_LINEAGE_INTERPRETATION.md](TARGET_LINEAGE_INTERPRETATION.md). Expression programs support lineage hypotheses; they do not disclose exact anonymous cell identities, establish a tissue of origin, or prove that a state transition happened.

## Reproduce from Command Prompt

Run from the project root with the existing competition environment and the unchanged input ZIP:

```cmd
.\.competition_env\Scripts\python.exe -m unittest test_target_archive -v
.\.competition_env\Scripts\python.exe analyze_target_archive.py --seed 20261004 --threads 12
.\.competition_env\Scripts\python.exe validate_target_archive.py
```

The validator does not call the analyzer's count or marker functions. It rereads the archive and original matrices, repeats all 16,560,000 exact count comparisons, checks every original cell ID/guide/depth, recomputes all 1,800 gene/group marker rows with sklearn ROC AUC, and checks every broader marker summary. It independently reranks 60 random cells plus all discordant cells, verifies saved centroid assignments and PCA projection, checks metrics/composition and earlier full-gene assignments, and verifies script/output hashes. Focused tests cover incomplete/duplicate target headers, ambiguous ZIP members, tied ranks, measured zeros, infinite/undefined ratios and full-library depth adjustment.

## Output files

- `archive_validation.json`: complete input schema/count audit and row-identity policy.
- `analysis_report.json`: parameters, convergence histories, versions, marker threshold counts and checksums.
- `validation.json`: independent verification evidence.
- `per_cell_target300_clusters.csv`: all 55,200 original IDs and assignments/coordinates/depth.
- `cells_disagreeing_with_context.csv`: all discordant cells and their earlier full-gene cluster.
- `cluster_composition.csv`, `cluster_metrics.csv`: complete contingency and K/seed comparisons.
- `target300_cluster_markers.csv`, `target300_context_markers.csv`: all 300 genes in each of three groups, with descriptive specificity metrics.
- `top10_enriched_cluster_markers.csv`: top qualifying markers per blind cluster.
- `target300_expression_summary.csv`: all 900 original-context/gene summaries.
- `broader_lineage_expression_summary.csv`: full-matrix lineage/proliferation checks.
- `marker_depth_quintiles.csv`: count-depth sensitivity summaries.
- `gene_axis.csv`, `pca_variance.csv`, `target300_rank_models.npz`: axes and frozen models.
- `target300_clusters_markers.png/.pdf`, `broader_lineage_markers.png/.pdf`: scientific figures.

No submission package is generated or uploaded by this analysis.
