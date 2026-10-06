# Predicting other genes from the 300-gene panel

This experiment fits three separate models: one for A, one for B and one for C. Each uses the measured counts of the 300 genes in `pert_counts.csv` to predict the other 18,233 genes. The fitted models and reconstructions are in `analysis/gene300_models/`; measured input counts are retained alongside every reconstructed cell.

These data are non-targeting control cells. The models learn associations between expression measurements in the same cell. Their control-cell accuracy does not establish that replacing one input value with a simulated CRISPRi knockdown predicts the biological response to that intervention.

## 1. Inputs and alignment

The original `context_A.h5ad`, `context_B.h5ad` and `context_C.h5ad` supply 18,400 cells each, with the same 18,533-gene axis. The existing `load_control` helper checks every stored count for finite, nonnegative integer values, canonicalizes the sparse matrices and verifies context labels and the 46 non-targeting guides with 400 cells each.

The 300 input columns follow `pert_counts.csv`; the output columns follow `gene_names.csv` with those 300 removed. Input genes never appear in output performance scores. The trainer rechecks all 16,560,000 values in the user-supplied `vcc_2026_300_target_matrix.zip` against the corresponding original H5AD measurements. The ZIP has no cell IDs, so source row order is essential. Original H5AD cell IDs are retained in split tables and prediction files. Input and output hashes are recorded in `run_manifest.json`.

## 2. Hold out entire guide groups

Within each context, seed 20261006 shuffles the sorted list of 46 control-guide IDs:

| Role | Guide groups | Cells | Use |
|---|---:|---:|---|
| Training | 32 | 12,800 | Fit features, expression coefficients and count calibration |
| Validation | 7 | 2,800 | Choose feature family and regularization |
| Test | 7 | 2,800 | Report final held-out performance |

A guide belongs to exactly one split. Neither cell IDs nor guide IDs are model features. There are no donor labels in these files, so this measures generalization to held-out control-guide groups within the same cell type, not unseen donors or cell types. `cell_splits.csv` lists every assignment.

All hyperparameters are written to `selected_A/B/C.json` before evaluating test expression. Feature means, scales and calibration parameters for the evaluation models use training rows only. Production models are then refitted using all 18,400 controls in their own context, retaining the selected hyperparameters. Production reconstructions are explicitly marked **in-sample**, and their accuracy is not reported as held-out accuracy.

## 3. Compare features

Each candidate uses only the 300 measured counts:

- **Rank:** average tied ranks within each cell, centered and divided by their Euclidean norm. Measured zeros are valid tied values. A constant input vector maps to a zero rank vector.
- **Magnitude:** `log1p` of the 300 raw counts.
- **Hybrid:** the 300 rank features and the 300 magnitude features together.

Every family also receives `log1p(sum of the 300 measured counts)` and the fraction of those 300 genes detected. Rank and magnitude models therefore have 302 features; hybrid has 602. Ranking and rank correlations add no pseudocount. `log1p` is a separate variance-reducing transform for regression; it is not used to alter values before computing rank correlations.

Feature centering and scaling use training cells only. Zero-variance features receive scale 1. The full-cell library size and expression of any of the 18,233 output genes are never supplied as inference inputs.

## 4. Expression regression

Training labels are `log1p(CPM)` for each output gene, where CPM = original count / original total counts across all 18,533 genes * 1,000,000. The full library total is used to construct supervised labels and evaluate predictions, not as an input feature.

Multi-output ridge regression minimizes `mean squared residual + alpha * squared coefficient norm`, with an unpenalized intercept and alpha in 0.01, 0.1, 1 and 10. A spectral solution fits all output genes without separate iterative optimizations. Validation loss is calculated from sufficient statistics; a unit test checks it against explicitly generated predictions.

Selection minimizes validation log-CPM MSE for genes with training mean CPM at least 5. This mask is based only on training expression and prevents never-expressed genes from dominating the comparison. All 18,233 output genes are still predicted and evaluated. Selection uses unconstrained regression outputs; reported test predictions clip log-CPM to the physical range 0 to `log1p(1,000,000)`.

## 5. Recover expected raw counts without knowing the full library

A separate ridge model, using the magnitude features of the same 300 genes, predicts `log1p(full library) - log1p(measured 300-gene total)`. Its regularization is selected on validation loss. Training residuals provide a multiplicative smearing factor when converting back from logarithms; training ratio bounds plus a one-unit margin limit extrapolation.

For each output gene, a smearing factor calculated from its training residuals corrects the inverse transformation of predicted log-CPM. Nonnegative corrected values define the relative expression profile across the 18,233 output genes. The predicted full library minus the observed 300-gene total is allocated across that profile. Consequently, predicted other-gene totals plus the original measured input counts equal the inferred full library, within float32 storage precision.

These outputs are **floating expected raw counts**, not sampled integer single-cell realizations. They are expression reconstruction artifacts, not `.vcc` competition submission files. The inferred library is an estimate from the 300 inputs, and its held-out relative error is reported separately.

## 6. Evaluate useful signal

`test_summary.csv` summarizes each context; `per_gene_test_metrics.csv` has all 54,699 context/output-gene combinations, and `per_cell_test_metrics.csv` contains the 8,400 held-out cells.

Compare against both context-average expression and a depth-only baseline. The depth-only raw-count predictor multiplies each test cell's measured 300-gene total by each output gene's training aggregate count divided by the training aggregate 300-gene count. It captures library-depth information without learning cell-state associations. The normalized-expression baseline is each gene's training mean log-CPM.

Skill is `1 - model MSE / baseline MSE`; positive values improve on the stated baseline, zero matches it and negative values underperform it. Per-gene normalized Spearman correlation ranks values across test cells with average ties; zeros remain valid. Constant vectors give undefined correlation, saved as blank CSV values. Per-cell profile correlations rank genes within each cell and are compared with the average profile because a high profile correlation alone can reflect common abundant genes.

The normalized correlations describe the predicted log-CPM layer before raw-count allocation. They are distinct from the raw-count Pearson correlation and raw-count reconstruction errors. See `GENE300_MODEL_RESULTS.md` for the actual results.

## 7. Saved artifacts

All paths below are relative to `analysis/gene300_models/`:

| File | Contents |
|---|---|
| `evaluation_model_A/B/C.npz` | Training-only models used for held-out evaluation |
| `model_A/B/C.npz` | Selected models refitted on all controls of their context |
| `heldout_predictions_A/B/C.h5ad` | 2,800 test cells x 18,233 expected-count outputs; predicted log-CPM layer |
| `reconstructed_A/B/C.h5ad` | All 18,400 cells x 18,233 outputs, marked in-sample |
| `candidate_scores.csv`, `depth_candidates.csv` | Every validation candidate |
| `cell_splits.csv` | Original IDs, guide groups, source rows and split assignments |
| `test_summary.csv` | Separate A/B/C held-out performance |
| `per_gene_test_metrics.csv`, `per_cell_test_metrics.csv` | Detailed held-out diagnostics |
| `run_manifest.json`, `validation.json` | Reproducibility and independent audit records |

H5AD `X` contains only the 18,233 predicted genes, named on the `var` index. `obsm['measured_300_raw_counts']` preserves the 300 measured inputs in `uns['input_genes']` order. `obs` preserves source cell IDs, contexts and guides, and adds prediction scope and inferred library size. The held-out files also contain `layers['predicted_log1p_CPM']`. Combining the output axis and named input axis recovers the full 18,533-gene panel without predicting the measured inputs.

Weights are NPZ files with numeric arrays and Unicode gene names, read with pickle disabled. Large model files and H5AD reconstructions stay local and are reproducible; small audit tables, code and documentation are committed.

## 8. Reproduce in Windows Command Prompt

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe test_gene300_models.py
.\.competition_env\Scripts\python.exe train_gene300_models.py --threads 8
.\.competition_env\Scripts\python.exe validate_gene300_models.py
.\.competition_env\Scripts\python.exe smoke_gene300_inference.py
.\.competition_env\Scripts\python.exe document_gene300_models.py
```

The existing `requirements-competition.txt` provides the installed dependency versions; this experiment requires no new packages. A full training run overwrites this experiment's generated artifacts. It does not modify any competition submission.

Apply a saved model to the supplied 300-gene ZIP:

```bat
.\.competition_env\Scripts\python.exe predict_gene300_models.py --context A --input vcc_2026_300_target_matrix.zip --output analysis\gene300_models\example_A.h5ad
```

Repeat with B or C to use the corresponding model. A CSV must contain the 300 named gene columns and may include `context`, `cell_id` and `ntc_id`. Column order is aligned by names. Missing genes, unexpected gene columns, duplicate names, negative/noninteger/nonfinite counts and duplicate cell IDs are rejected. Output files are not overwritten by the inference CLI. ZIP inputs without cell IDs receive explicit input-row IDs.

## 9. Independent verification

`test_gene300_models.py` checks tied ranks, constant inputs, deterministic guide isolation, the ridge solution against augmented least squares, validation sufficient statistics, held-out label/depth isolation, serialization, batch invariance, count conservation and named input alignment.

`validate_gene300_models.py` checks all source/output hashes and source cell/gene identities, every split, and saved feature statistics. It independently solves the normal equations for 24 genes per evaluation and production model, and checks all depth coefficients. It recalculates every held-out gene's raw and log-CPM MSE and both baseline errors, replays all 8,400 held-out cells across all output genes, scans every full reconstruction value for finite/nonnegative counts and count conservation, and replays 24 full-data rows per context. Its audit record states exactly which checks are exhaustive and which are sampled.

`smoke_gene300_inference.py` additionally runs the standalone prediction CLI on real measured cells from A/B/C with the input columns deliberately reversed. It verifies saved counts, gene order and preserved cell IDs; results are in `inference_smoke.json`.
