# Neural networks for the 300-gene prediction task

This experiment compares four neural configurations separately in A, B and C. Each receives only the 300 measured input genes and predicts the remaining 18,233 genes. Each configuration is trained twice: once for normalized expression and once directly for raw counts, giving 24 candidate networks. The best validation configuration for each context/objective is saved and refitted on all controls.

## 1. Same measured cells and gene axes

Training uses the original `context_A/B/C.h5ad` control matrices: 18,400 cells and 18,533 genes in each context. The existing verified loader checks all raw counts, context labels, gene order and guide balance. The input axis follows the 300 symbols in `pert_counts.csv`; the output axis is `gene_names.csv` with those inputs removed. Predictions never include input genes in their error scores.

This is observational reconstruction of control-cell expression. It does not provide experimentally validated responses to CRISPRi interventions.

## 2. Same guide-group comparison with ridge

The seed and split are the same as the ridge study: 32 control-guide groups/12,800 cells for training, seven/2,800 for validation and seven/2,800 for test. The neural trainer explicitly checks its split against the saved ridge split table. Feature statistics, baseline profiles and training labels use training rows only; validation chooses architecture and epoch.

The test set was already reported in the ridge study. The neural comparison is therefore a comparative evaluation on that existing test set, not a new independent test. Neural weights, architecture choices and stopping epochs are locked before scoring these test cells. All four frozen configurations are then evaluated once; test scores do not change the selection.

There are no donor labels, so held-out guide-group performance does not establish generalization to unseen donors or new cell types.

## 3. Input features and neural structures

The magnitude input consists of `log1p` of the 300 observed counts, `log1p` of their total and their detection fraction (302 features). The hybrid input adds 300 within-cell rank features (602 features). Rank features retain zeros, use average ties, center the ranks and divide by their Euclidean norm. No pseudocount is added for ranks or rank correlations. Constant rank vectors map to zeros. Feature means and scales are fitted on training rows only.

| Configuration | Learned network | Activation | Dropout | Input family |
|---|---|---|---:|---|
| Shallow | 302 -> 128 -> 18,233 | ReLU | 0.10 | Magnitude |
| Deep | 302 -> 256 -> 128 -> 18,233 | SiLU | 0.10 | Magnitude |
| Bottleneck | 302 -> 128 -> 32 -> 18,233 | SiLU | 0.05 | Magnitude |
| Residual | 602 -> 256 -> 128 -> 18,233 plus a linear 602 -> 32 -> 18,233 branch | GELU | 0.10 | Hybrid |

For normalized expression, the residual configuration additionally starts from the training-only ridge prediction and learns a correction. Other normalized configurations start from each gene's training mean. All raw configurations start from the same depth-only raw-count baseline. The residual configuration thus differs in input features, the learned branch and the normalized baseline; this is a comparison of useful configurations, not an isolated causal test of architecture depth.

Output weights and biases initially equal zero, so epoch zero exactly recovers each configuration's baseline. Epoch zero is included in validation selection; a configuration that fails to improve its baseline retains that baseline. Hidden layers start from PyTorch's seeded initialization.

## 4. Two distinct prediction objectives

**Normalized:** targets are `log1p(CPM)`, using original full-library totals to form supervised labels. Training minimizes MSE across all outputs. Validation selects MSE on genes with training mean CPM at least 5. Inference clips predictions to 0 through `log1p(1,000,000)`.

**Raw:** the fixed baseline for output gene j is the cell's measured 300-gene total multiplied by the aggregate training counts for gene j divided by the aggregate training 300-gene counts. The network learns a gene-specific log-ratio correction, and the expected count is `baseline * exp(clamp(correction, -3, 3))`. Training and validation compare these expected counts directly with original raw counts across all 18,233 outputs. The loss is divided by a single training-label mean-square scale for stable optimization; this does not change the objective's relative weighting of genes.

The raw head's outputs are nonnegative floating expected counts. Genes with zero counts throughout training retain a zero baseline. An all-zero measured panel also produces zero raw outputs, a boundary case outside the observed control-panel range. The normalized head is saved separately and is not converted into raw counts using the earlier smearing procedure.

Full-library totals and any of the other 18,233 genes are never provided as inference inputs. Raw predicted library size is the sum of the inferred other-gene counts plus the original measured 300-gene counts.

## 5. Optimization and early stopping

The prespecified search uses up to 40 epochs, batch size 256, AdamW with learning rate 0.001 and weight decay 0.0001, gradient norm capped at 5, and seven epochs of patience. Validation must improve by more than 0.000001 to replace the saved best epoch. Every epoch's train and validation errors are recorded. The architecture with lowest best validation MSE wins separately for each context and objective.

GPU training uses float16 autocast for network operations, float32 losses and outputs, and gradient scaling. Validation and saved inference use float32 network computation. Seeds are recorded for every candidate. This is a first comparison with one seed per configuration; small differences between configurations should be checked with repeated seeds in a subsequent study.

Selected production models are initialized with the selected seed and refitted on all 18,400 controls for the selected number of epochs, without using test error to stop training. These production reconstructions are in-sample and are labelled accordingly.

## 6. Saved models and files

All experiment artifacts are under `analysis/gene300_neural/`:

| File | Meaning |
|---|---|
| `candidate_A/B/C_normalized/raw_ARCH.npz` | All 24 best-epoch candidate checkpoints |
| `candidate_scores.csv` | Validation scores, epochs, seeds, parameter counts and model names |
| `training_history.csv` | Every recorded training/validation epoch |
| `selected_A/B/C.json` | Settings locked before the neural test evaluation |
| `evaluation_normalized/raw_A/B/C.npz` | Selected training-only models |
| `model_normalized/raw_A/B/C.npz` | Selected models refitted on all controls |
| `test_architecture_scores.csv` | One held-out evaluation of every frozen candidate |
| `heldout_normalized_A/B/C.h5ad` | Selected normalized predictions for 2,800 held-out cells each |
| `heldout_raw_A/B/C.h5ad` | Selected expected raw counts for the same test cells |
| `reconstructed_raw_A/B/C.h5ad` | Expected raw counts for all 18,400 cells in each context |
| `per_gene_test_metrics.csv` | 109,398 selected-model context/objective/gene rows |
| `per_cell_test_metrics.csv` | 16,800 selected-model context/objective/cell rows |
| `run_manifest.json`, `validation.json` | Versions, source/output hashes and independent audit coverage |

NPZ checkpoints store neural weights, preprocessing statistics, the fitted baseline, gene axes and scope as numeric/Unicode arrays. They are loaded with pickle disabled. Residual normalized checkpoints include their fixed ridge arrays and therefore do not need another weight file during inference. PyTorch is required to evaluate these neural weights.

H5AD `X` contains only the 18,233 output genes. **Normalized files contain log1p(CPM); raw files contain expected raw counts.** Units and objective are explicitly saved in `uns`. The measured input counts are retained in `obsm['measured_300_raw_counts']`, with names in `uns['input_genes']`. Original source cell IDs, context and guide metadata are preserved in `obs`. Full-data normalized reconstructions can be generated with the inference CLI when needed; the default experiment saves full-data raw reconstructions and both kinds of held-out predictions.

Large checkpoints/H5AD files remain local and reproducible. Code, Markdown, figures, candidate tables, test errors and audit records are versioned.

## 7. GPU runtime and Windows shared-drive setup

The project uses Python 3.14 and the official pinned CUDA wheel `torch==2.10.0+cu126`, installed into `.competition_env`. The RTX 2070 SUPER passed a real CUDA forward/backward calculation. [Official PyTorch installation instructions](https://pytorch.org/get-started/previous-versions/) list the wheel indexes; [NVIDIA's minor-version compatibility documentation](https://docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html) explains the limited compatibility of CUDA applications within one major driver family. The actual kernel check is the evidence for this machine, not merely its displayed CUDA version.

A Windows DLL-loading error occurred during shared-drive runtime initialization. `cache_neural_runtime.py` provides a local temporary copy of the installed package, excluding C++ headers/share files and bytecode. Every copied DLL is checksum-verified against the installed official wheel. No binaries, system drivers or paging settings are modified. `neural_runtime.py` activates that cache from a local experiment configuration; this machine-specific pointer is ignored by Git. The cache audit is retained. If the temporary cache is cleared, run the cache script again.

## 8. Reproduce in Windows Command Prompt

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m pip install torch==2.10.0 --index-url https://download.pytorch.org/whl/cu126
.\.competition_env\Scripts\python.exe cache_neural_runtime.py
.\.competition_env\Scripts\python.exe test_gene300_neural.py
.\.competition_env\Scripts\python.exe train_gene300_neural.py --device cuda --epochs 40 --patience 7 --batch-size 256
.\.competition_env\Scripts\python.exe validate_gene300_neural.py
.\.competition_env\Scripts\python.exe smoke_gene300_neural.py
.\.competition_env\Scripts\python.exe document_gene300_neural.py
```

`requirements-neural.txt` pins the complete dependency set. The comparison expects the previously generated ridge model weights and split table in `analysis/gene300_models/`; reproduce those first using `GENE300_MODELS.md` if absent. Training overwrites this neural experiment's generated outputs. It does not modify the ridge experiment or a competition submission.

Predict expected raw counts from the provided ZIP:

```bat
.\.competition_env\Scripts\python.exe predict_gene300_neural.py --context A --objective raw --input vcc_2026_300_target_matrix.zip --output analysis\gene300_neural\new_raw_A.h5ad
```

Use `--objective normalized` for log-CPM predictions, and B/C for their own context-specific models. The CLI accepts the same named-column CSV/ZIP schema as the ridge predictor, aligns gene columns by name, checks count validity and refuses to overwrite output files. CPU inference is the default; `--device cuda` enables GPU inference.

## 9. Verification

Five focused tests check all architectures' output dimensions and gradients, baseline initialization, raw nonnegativity, training-only preprocessing/baselines, all eight architecture/objective combinations against an independent NumPy forward calculation, NPZ replay, correct best-epoch selection, and test-label isolation during training.

The independent validator checks all selected artifact/source hashes, exact guide splits and fitted feature statistics, validation-only architecture/epoch selection, preserved input/output gene axes and source cell IDs, training-only raw baselines and fixed ridge arrays, and explicit NumPy neural inference. It recomputes every selected model's gene and cell MSE and corresponding baseline errors from original counts. Every saved expression value is scanned for finite/nonnegative values. Selected evaluation/production models and full reconstructed files are independently replayed for 32 rows each. Its JSON states the exhaustive and sampled coverage.

`smoke_gene300_neural.py` also checks all six standalone context/objective commands on real measured counts with reversed input-column order, preserving original cell IDs and output gene names. Its results are saved in `inference_smoke.json`.
