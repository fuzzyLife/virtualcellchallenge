# External controls mapped to a competition-only rank space

This analysis places external K562, HCT116, HEK293T, H1 and CD4 control observations relative to supplied A/B/C controls. It preserves the raw-count, average-tied-rank correlation strategy and records both external unit-level dispersion and pooled-profile similarity.

## Why the gene axis changes

The previous competition-only clustering used every one of the 18,533 measured genes. None of the external sources measures all of them:

| Reference | Competition genes measured |
| --- | ---: |
| K562 | 7,681 |
| HCT116 | 18,106 |
| HEK293T | 18,106 |
| H1 | 18,077 |
| CD4 | 17,771 |

An exact projection into the previous complete-data rank coordinates cannot be calculated for a reference with missing genes. Computing ranks on a shorter gene axis also changes their values. This study therefore uses the **same rank-space construction on the 7,447 jointly measured genes**, fits it on all 55,200 competition cells, and freezes it before projecting the references. Its PCA coordinates differ from the previous 18,533-gene plot.

Symbols are matched exactly, retaining competition gene order. Duplicate measured source symbols are summed. Every common-axis gene must be measured in every source. An omitted entry in a sparse count row is a genuine measured zero; a symbol absent from the source's feature metadata is unmeasured and excluded from the common axis. No missing symbol is filled with zero. Every ranked value must be finite; constant vectors and missing values stop the analysis.

## 1. Inspect and recover the existing controls

Root scripts `prepare_reference_controls.py`, `analyze_cell_reference_correlations.py`, `cluster_full_rank_matrix.py` and `extract_h1_responses.py` were read before implementation. The same source releases and selected control rows as the earlier reference study are reused:

| Reference | Observations mapped | Observation type | Represented cells |
| --- | ---: | --- | ---: |
| K562 | 503 | Reconstructed non-targeting pseudobulk count sums | 75,094 |
| HCT116 | 891 | Passing non-targeting single cells in Batch1 | 891 |
| HEK293T | 1,079 | Passing non-targeting single cells in Batch1 | 1,079 |
| H1 | 1,536 | 32 source-ordered non-targeting cells from each of 48 batches | 1,536 |
| CD4 | 384 | Raw non-targeting pseudobulk sums, 32 per donor/condition/run stratum | 34,620 |

There are **4,393 external observations**. The aggregate source releases do not provide each represented K562/CD4 cell separately; their rows must not be described as individual cells. CD4 Rest, Stim8hr and Stim48hr conditions are retained individually.

`extract_external_mapping_controls.py`:

1. Checks the eight prepared pooled-profile hashes and the five selected-row hashes against the verified earlier reference analysis.
2. Reads K562 original raw mean-expression vectors and multiplies each by its filtered cell membership; collapses duplicate symbols.
3. Reads the two pinned Orion parquet files and original gene-token metadata; extracts the exact passing selected control cells. Sparse token omissions represent measured zero counts.
4. Reads the SHA-verified H1 selected-count CSR cache, selects only its 1,536 non-targeting rows, and verifies the original source-row identifiers.
5. Opens CD4 through `CachedRangeReader`, an offline HDF5 file interface. Each consumed 256-KiB block must match its cached SHA-256 sidecar, and the cached file size/ETag/URL must match the verified source identity. Missing or corrupt blocks stop extraction. No remote request is made; consumed bytes remain bounded to 2 GiB.
6. Checks that the sum of reconstructed external observations recovers each previously verified pooled reference on the common axis. Saves the aligned sparse counts, source metadata, ordered gene axis and detailed extraction provenance under ignored `reference/external_mapping/`.

Data sources and their original selection limitations are detailed in [CELL_REFERENCE_CORRELATIONS.md](CELL_REFERENCE_CORRELATIONS.md). This maps the previously selected reference subsets, not the entire public atlases.

## 2. Fit on competition controls and freeze

`map_external_rank_controls.py` verifies original H5AD hashes against the completed full-rank clustering report. It keeps all 18,400 original control cells in each context and reranks their raw counts on the common axis:

```
z = (average_rank(raw_counts) - mean_rank) / norm(centered_ranks)
Spearman(cell_i, cell_j) = dot(z_i, z_j)
```

There is no pseudocount, log, CPM transform, gene variance scaling or missing-value imputation. Measured zeros remain valid ties. Ranks are calculated in double precision and stored in float32, as in the completed full-gene study.

The existing spherical clustering algorithm is reused: label-blind MiniBatchKMeans initialization, followed by full-cell/full-gene spherical refinement until assignments cease to change. K=3 is fit with seeds 20261004 and 20261005, and agreement is reported. The A/B/C names are attached to fitted clusters afterward through the one-to-one assignment maximizing the context-versus-cluster count table. This naming does not alter the fitted centroids.

Eight seeded randomized PCs are fitted solely to the competition rank vectors; coordinates use explicit projection onto the saved components. Both centroids and PCA are frozen before external observations enter the mapping. External expression does not influence clustering or PCA. The common gene axis depends on source feature coverage rather than external expression values.

## 3. Project and quantify each external observation

Each external observation is reranked on the same complete common axis. Its dot products with the three frozen unit centroids determine its closest A/B/C context. The margin is the largest similarity minus the second largest. These are similarities to mean-rank centroids, not correlations to raw cells that have been reranked, and are not calibrated probabilities.

The mean of the actual Spearman correlations to every competition cell in a context is calculated exactly up to float32 storage precision by dotting the external unit-rank vector with that context's mean unit-rank vector. All 18,400 peers per context are included, without materializing a 4,393-by-55,200 pairwise matrix.

The external PCA coordinates are `(z_external - competition_PCA_mean) @ competition_components.T`; there is no joint PCA refit. Two PCs are plotted for readability, with explained-variance fractions. Full-gene centroid assignment uses all 7,447 features, not the plotted PCs.

Two distribution-shift diagnostics are retained:

- Compare each observation's best centroid similarity with the assigned context's fifth percentile among the original competition controls.
- Compute squared Mahalanobis distance in all eight PCs using the original context's PC mean and empirical covariance pseudoinverse; compare with its original empirical 95th percentile. Distances to all three contexts are also retained, with a flag for observations outside every context's 95th-percentile region.

These thresholds describe the supplied control clouds. They are not calibrated tests of biological identity, especially for pseudobulks whose depth and pooling differ from individual competition cells. External observations are displayed at their actual projected positions even when outside the competition clouds.

## 4. Check pooled-profile sensitivity

The five existing pooled references and three CD4 condition-specific pools are mapped through the same frozen model. The three competition raw-count pools are also projected as visible plus-marker anchors.

Separately calculate ordinary Spearman correlations between every raw pooled profile and the three raw competition pools, with the same 7,447 measured genes and tied zeros. Report `nearest_pooled_context` alongside the unit-centroid-based `nearest_context`. Pooling counts before ranking differs from averaging cell rank vectors; the two comparisons may disagree. No conclusion is selected by comparing unaligned gene axes or by treating pooled observations as individual cells.

Assignment fractions count source observations, not represented cells. Per-observation membership, source row, guide, batch, culture condition, donor and run metadata remain in the CSV, as available. CD4 conditions have separate summaries.

## 5. Reproduce in Command Prompt

From the project root with the existing source caches:

```cmd
.\.competition_env\Scripts\python.exe -m unittest test_external_rank_mapping -v
.\.competition_env\Scripts\python.exe extract_external_mapping_controls.py
.\.competition_env\Scripts\python.exe map_external_rank_controls.py --seed 20261004 --threads 12
.\.competition_env\Scripts\python.exe validate_external_rank_mapping.py
```

On a fresh checkout, first restore the original competition matrices and follow [CELL_REFERENCE_CORRELATIONS.md](CELL_REFERENCE_CORRELATIONS.md) to download and prepare the pinned public controls, and `extract_h1_responses.py` to recreate the H1 selected-count cache. The previous correlation and full-rank reports must be present for the provenance checks. Install the pinned libraries in `requirements-competition.txt`. The mapping run itself requires no additional download. Original data and prediction files are not edited.

All scripts stay in the root. The extraction cache is ignored by Git; the frozen model, per-observation coordinates, gene axis, figures, provenance and summaries are tracked. Reported settings, software versions, hashes and convergence histories permit reproduction. Budget at least 12 GiB of free memory; the observed workstation has about 128 GiB RAM.

## Outputs

Tracked files in `analysis/external_rank_mapping/`:

| File | Contents |
| --- | --- |
| `gene_axis.csv` | Ordered 7,447 jointly measured gene symbols |
| `competition_cells_shared_space.csv` | All original 55,200 IDs, context/guide metadata, common-axis count depth/detection, fitted cluster and context, three centroid similarities and eight projected PCs |
| `external_control_units_mapped.csv` | All 4,393 external units, original source metadata, closest context, three similarities, three exact mean cell-level Spearman correlations, margins, PCs and shift diagnostics |
| `external_mapping_summary.csv` | Observation-level A/B/C fractions, represented membership, median similarities/margins and shift fractions for every source and CD4 condition |
| `pooled_profiles_mapped.csv` | Eight external and three competition pools, frozen coordinates, centroid assignment and independent raw pooled-profile Spearman comparison |
| `pooled_raw_counts.csv` | Actual common-axis raw-count vectors for the 11 pooled profiles |
| `competition_cluster_composition.csv` | Blind shared-axis cluster composition, evaluated afterward against A/B/C |
| `pca_variance.csv` | Explained variance of eight competition-only PCs |
| `frozen_shared_rank_model.npz` | Unit centroids, competition mean-rank vectors, PCA parameters and context-specific shift calibration; pickle-free |
| `external_control_extraction.json` | Source/control selection, count-cache hashes and every consumed CD4 block checksum |
| `analysis_report.json` | Fit scope, dimensions, seeds, histories, source coverage, versions, script and output hashes |
| `external_rank_projection.png`, `.pdf` | All competition cells and external observations, with separate source panels and pooled anchors |
| `external_mapping_summary.png`, `.pdf` | Pooled-profile Spearman heatmap and observation-level assignment fractions |
| `validation.json` | Independent verification results |

[EXTERNAL_RANK_MAPPING_RESULTS.md](EXTERNAL_RANK_MAPPING_RESULTS.md) records the measured results.

## Verification

Focused tests cover summing duplicate gene symbols, retaining measured zero entries, rejecting missing gene axes/nonfinite observations, checking frozen projection against SciPy Spearman, verifying that row order cannot refit the model, and rejecting missing/corrupt offline source blocks.

The independent validator checks every competition ID and external source-row ID, source/output/script hashes, the measured common axis, all assignment-fraction tables and cluster agreement. It reranks 30 actual competition cells and 25 external observations reconstructed independently from original source rows, checks frozen centroid assignments and PCA coordinates, and verifies mean correlations and eight-PC distances. It independently checks all 11 pooled profiles and 33 actual pooled Spearman correlations. The extraction also proves that source-unit sums reproduce previously verified pools.

Method reference: [PCA transform documentation](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html). Context proximity characterizes basal expression; it does not establish cell identities or transferable CRISPRi responses.
