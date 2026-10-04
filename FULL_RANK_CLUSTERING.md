# Full-matrix Spearman clustering

This analysis asks whether the supplied A/B/C non-targeting controls form separate clusters when all cells and all measured genes are considered together. The root scripts are `cluster_full_rank_matrix.py`, `test_full_rank_clustering.py`, and `validate_full_rank_clustering.py`.

## Inputs and scope

- `context_A.h5ad`, `context_B.h5ad`, `context_C.h5ad`: 18,400 original controls per context, 55,200 cells total.
- `gene_names.csv`: the complete ordered axis of 18,533 genes. No restriction to the 300 competition targets or to external-reference gene coverage.
- Original cell IDs, context labels, `target_gene`, and `ntc_id` are retained. There are 46 non-targeting guide IDs per context.
- Input SHA-256 values must match the previously verified target-rank study. Each cell-by-gene count is checked for finiteness, nonnegativity, and integrality. All source gene orders must match.

Only these supplied controls are clustered. Neither external cell lines nor predicted perturbations are included.

## Rank geometry

For each cell, rank its raw counts across all 18,533 genes using average ranks for ties. Measured zeros remain valid tied values. Subtract the row's mean rank and divide by its Euclidean norm:

```
z_i = (rank(x_i) - mean(rank(x_i))) / norm(rank(x_i) - mean(rank(x_i)))
rho(i,j) = dot(z_i, z_j)
Spearman distance(i,j) = 1 - rho(i,j)
```

There is no added constant, log transform, CPM normalization, selection of variable genes, imputation, or gene standardization. Unit scaling of centered ranks is the algebra needed to calculate correlation, rather than a count transformation. The ranks are calculated in double precision and stored in float32; reported full-matrix coefficients therefore have floating-point error, checked against double-precision SciPy results.

The earlier pairwise-valid rule is preserved: every input here is finite, so every cell pair has the same 18,533 valid gene pairs. If any missing or infinite value appears, this whole-matrix implementation stops. Pair-specific reranking would require a different representation; missing entries are never replaced with zeros. Constant cell vectors also stop the analysis because their correlations are undefined.

## Clustering without context labels

1. Concatenate the three original matrices in A, B, C order, retaining every cell. Construct the full 55,200-by-18,533 unit-rank matrix in memory (about 4.1 GB).
2. For K=2,3,4,5,6, initialize from scikit-learn MiniBatchKMeans on the full gene matrix (`n_init=3`, batch size 1,024, initialization sample 8,192, maximum 30 passes, patience 15). A fixed seed makes sampling reproducible.
3. Refine on **all cells and all genes** with spherical k-means. Assign each cell to the unit centroid with the greatest rank-feature dot product; sum assigned cell vectors and normalize the centroid. Stop when assignments cease to change; reject a fit that does not converge within 100 iterations or creates an empty cluster. This maximizes the sum of cell-to-unit-centroid correlations. A centroid is a mean of cell rank vectors, not a raw-count cell that is reranked.
4. Repeat K=3 with an independent seed. Compare the two assignments using adjusted Rand index (ARI).
5. Calculate silhouettes using exact `1 - rho` distances for the same uniformly sampled 1,500 cells in every fit. Sampling limits the quadratic evaluation cost; **clustering itself is not sampled or projected**. Evaluate the supplied context labels on the same sample for comparison. Report full cluster sizes and the number represented in the silhouette sample; tiny clusters may be omitted from that sample, which limits the comparison. Exact score ties are resolved toward smaller K. Report all candidate K values; the largest sampled silhouette in this finite range is a diagnostic, not proof of a biological cluster count.
6. Only after fitting, compare clusters with the A/B/C labels: complete count and proportion tables, ARI, normalized mutual information (NMI), and cluster purity. A purity of one indicates every cluster contains a single context, but purity tends to increase as more clusters are requested.

The program never passes context or guide labels to the fitter. Arbitrary cluster numbers carry no biological meaning. The initializer is stochastic and the refinement finds a local optimum; the repeat checks sensitivity to initialization but does not exhaust all local optima.

## Dispersion, confounding and visualization

- Every cell has its assigned-centroid correlation and margin over the next centroid. Small margins identify ambiguous assignments.
- Each cell's mean correlation to every A/B/C peer is computed from the full gene ranks and context sums, excluding self. This is an exact all-peer calculation up to stored float32 precision; no dense 55,200-square correlation matrix is needed.
- Eight PCs are fit by seeded randomized PCA to the full rank matrix, solely for display and diagnostics. Neither the PCA scores nor the two-dimensional plot define the clusters. Projection overlap cannot override full-matrix assignments; PC1/PC2 variance fractions are shown.
- Within each context, relate PCs to raw count depth and detected-gene count using Spearman correlation. Tabulate cluster membership within count-depth quintiles and within every non-targeting guide; report guide/cluster NMI within each context.
- A difference between contexts may reflect biology, assay or batch differences. These diagnostics do not identify the anonymous contexts or demonstrate successful perturbation transfer. Rank correlation still depends on dropout and tied zeros associated with count depth.

## Reproduce in Command Prompt

From the project root, using the existing environment:

```cmd
.\.competition_env\Scripts\python.exe -m pip install -r requirements-competition.txt
.\.competition_env\Scripts\python.exe -m unittest test_full_rank_clustering -v
.\.competition_env\Scripts\python.exe cluster_full_rank_matrix.py --seed 20261004 --threads 12
.\.competition_env\Scripts\python.exe validate_full_rank_clustering.py
```

The installed versions, input and output SHA-256 values, algorithm settings, seeds, convergence histories, and elapsed time are recorded in `analysis/full_rank_clustering/analysis_report.json`. The analyzed workstation has about 128 GiB RAM. Budget at least 16 GiB of available memory for the full matrix, PCA copies and temporary batches; runtime depends on CPU, BLAS and storage. The script overwrites only this analysis folder and its result Markdown when rerun.

## Outputs

All numeric artifacts are under `analysis/full_rank_clustering/`:

| File | Contents |
| --- | --- |
| `per_cell_full_rank_clusters.csv` | One row per original cell; original identifiers, raw count depth, detected genes, zero fraction, labels for all fits, K=3 centroid correlation/margin, exact means to all A/B/C peers and eight PCs |
| `gene_axis.csv` | All 18,533 gene symbols in source order |
| `cluster_composition_k3.csv` | Full A/B/C versus cluster count table and fractions within context/cluster |
| `cluster_composition_all_fits.csv` | Full context-versus-cluster counts for every fit, including singleton clusters |
| `within_context_dispersion.csv` | Context-wise means, SDs and quantiles of depth, detection and K=3 correlation/margin |
| `cluster_metrics.csv` | K sensitivity, repeat seed, convergence, sampled silhouettes, ARI/NMI/purity and fit objective |
| `silhouette_sample.csv` | The fixed 1,500 sampled original IDs and concatenated row indices |
| `all_peer_context_mean_correlations.csv` | Full-gene A/B/C mean pairwise Spearman correlations, excluding self |
| `depth_quintile_composition.csv` | Cluster membership in each context's five count-depth strata |
| `guide_composition.csv` | Membership for every guide in each context |
| `pc_depth_correlations.csv` | Within-context rank correlations of each PC with count depth and gene detection |
| `pca_variance.csv` | Variance explained by the eight display PCs |
| `cluster_models.npz` | Six fits' unit centroids and PCA mean/components; no pickle objects |
| `full_rank_clusters.png`, `.pdf` | All-cell PCA plots colored by context, blind cluster and depth, plus composition and sensitivity diagnostics |
| `analysis_report.json` | Provenance, software versions, hashes and refinement history |
| `validation.json` | Independent validation results |

`FULL_RANK_CLUSTERING_RESULTS.md` records the results. The original competition and submission files are not edited.

## Verification

The focused tests compare unit-rank dot products with SciPy Spearman for tied zeros, reject missing/constant inputs, and check that blind spherical refinement recovers known rank patterns with nondecreasing objective and consistent unit centroids.

The independent validator verifies every original cell ID, the gene axis, all file/script/input hashes, cluster counts and proportions, ARI/NMI/purity, convergence, and repeat agreement. It independently reranks 60 original cells in double precision and checks all six saved fits (360 assignments), original count-depth/detection summaries, centroid correlations and PCA coordinates. It also compares 105 actual-cell pairs with SciPy Spearman, independently reranks the complete 1,500-cell silhouette sample from the original H5AD files, and recomputes all six silhouette scores in double precision. Model correctness and context separation do not establish prediction accuracy on the leaderboard.

Method references: [MiniBatchKMeans documentation](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.MiniBatchKMeans.html), [silhouette documentation](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.silhouette_score.html).
