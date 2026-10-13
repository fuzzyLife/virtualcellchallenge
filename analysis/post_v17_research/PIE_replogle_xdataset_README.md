---
library_name: pytorch
license: other
license_name: pie-model-license
license_link: https://huggingface.co/arcinstitute/PIE_replogle_xdataset/blob/main/MODEL_LICENSE.md
tags:
  - biology
  - gene-expression
  - perturbation-prediction
  - perceiver
  - cross-dataset-generalization
datasets:
  - arcinstitute/PIE_tahoe100m
  - arcinstitute/PIE_jiang
  - arcinstitute/PIE_arc_vcc_25
  - arcinstitute/PIE_x_atlas_orion
  - arcinstitute/PIE_sources
---

# PIE: cross-dataset generalization to Replogle

PIE predicts the transcriptional response to a perturbation in a cell context using biological
knowledge embeddings, context control expression, and pooled response evidence from the training
set.

Code and CLI documentation: [ArcInstitute/PIE](https://github.com/ArcInstitute/pie)\
Datasets and knowledge sources: [PIE collection](https://huggingface.co/collections/arcinstitute/pie)

For every measured gene, it returns:

| Output | Meaning |
|---|---|
| `p_de` | Predicted probability of differential expression |
| `lfc_pred` | Predicted log2 fold change |
| `delta_p_pred` | Predicted change in mean expression on the preprocessed expression scale |

This repository contains one model from the `replogle_xdataset` experiment. It trains on Tahoe,
Jiang, ARC VCC 25, and Orion, and evaluates on Replogle without training or fine-tuning on
Replogle response labels. The benchmark measures transfer across datasets for perturbations
that are either seen or unseen in training.

## Checkpoint and files

```text
README.md
best_auprc.ckpt
config.yaml
data_stats.json
eval/
  test_seen/{metrics_best_auprc.csv,granular_best_auprc.csv}
  test_unseen/{metrics_best_auprc.csv,granular_best_auprc.csv}
```

`best_auprc.ckpt` is the Lightning checkpoint with the highest validation `binary_auprc`.
Validation uses the four training datasets. The checkpoint embeds its model configuration, pinned asset
references, and fitted data statistics; the adjacent YAML and JSON files make those settings
inspectable.

## Data and training

The saved run uses the following preprocessed datasets.

| Dataset | Role | Pinned HF revision |
|---|---|---|
| [Tahoe](https://huggingface.co/datasets/arcinstitute/PIE_tahoe100m) | Training and validation | `597d6ded9f9a148331388484ffb193489121c702` |
| [Jiang](https://huggingface.co/datasets/arcinstitute/PIE_jiang) | Training and validation | `99242356dfcf58e049d867703c5743df8d1790bf` |
| [ARC VCC 25](https://huggingface.co/datasets/arcinstitute/PIE_arc_vcc_25) | Training and validation | `e7dda6065316959958b87583ddf0a15f6e69c37d` |
| [Orion](https://huggingface.co/datasets/arcinstitute/PIE_x_atlas_orion) | Training and validation | `0950c8aa5ce5b33dc9d6f0c5eade3ae445ca0ffb` |
| [Replogle-Nadig essential](https://huggingface.co/datasets/arcinstitute/PIE_replogle_nadig_essential) | Evaluation; anchors gene-axis ordering | `20c9faef76fc96fdc809871340bd93499413dfe7` |

Replogle remains in the data configuration for evaluation, with no training or validation
pairs. Its response labels do not contribute to training or pooled response
evidence. Its gene axis anchors the union's ordering; Replogle evaluation reports predictions
for its 6,642 measured genes. Target-context control expression is an inference input, so
zero-shot here means no target response-label training or fine-tuning.

Canonical splits live in
[PIE splits](https://huggingface.co/datasets/arcinstitute/PIE_splits) under `replogle_xdataset/`.
The training split contains 64,260 context–perturbation pairs, and validation contains 957.

Knowledge sources are `esm2`, `ncbi_text`, `string_space`, `depmap_gene_effect`, `context_text`,
`perturbation_text`, `smiles`, `l1000_tas`, `prism_secondary`, and `jump_morphology`, plus
`gene_text` for gene queries. All come from
[PIE sources](https://huggingface.co/datasets/arcinstitute/PIE_sources) at revision
`cb1aaa4e7655605bdc70a9bd77bbd62016b8c7d7`.

## Evaluation

Both settings evaluate `best_auprc.ckpt` on HepG2, Jurkat, K562, and RPE1 in Replogle, with no
Replogle response-label training or fine-tuning. Seen/unseen refers to whether the perturbation
identifier occurs in the cross-dataset `train.json`.

| Setting | Target dataset | Perturbations | Test pairs |
|---|---|---|---:|
| `test_seen` | Replogle | Seen in training | 2,286 |
| `test_unseen` | Replogle | Unseen in training | 3,173 |

Values below are **mean ± one sample standard deviation across the four Replogle cell lines**
(`ddof=1`), with equal weight per cell line.

| Setting | AUPRC | Jaccard | Direction match | Spearman DE count | Spearman LFC | L1 discrimination |
|---|---|---|---|---|---|---|
| `test_seen` | 0.3802 ± 0.0242 | 0.1523 ± 0.0176 | 0.7401 ± 0.0236 | 0.6231 ± 0.0542 | 0.5130 ± 0.0421 | 0.6312 ± 0.0414 |
| `test_unseen` | 0.2463 ± 0.0560 | 0.1196 ± 0.0010 | 0.6466 ± 0.0121 | 0.5016 ± 0.0114 | 0.2698 ± 0.0614 | 0.5309 ± 0.0140 |

Metric names in the CSVs are `binary_auprc`, `sig_jaccard`, `direction_match`, `spearman_nsig`,
`spearman_lfc`, and `discrimination_score_l1`. Higher is better for all six.

## Download and use

Install PIE with `pip install arc-pie` (or, from a checkout of the PIE code repository,
`uv sync --frozen` for the pinned environment). Set the writable data, runs and cache
directories in `common.sh` (copy `common.sh.example` from the code repository):

```bash
pip install arc-pie
hf auth login

hf download arcinstitute/PIE_replogle_xdataset \
  best_auprc.ckpt config.yaml data_stats.json \
  --local-dir models/replogle_xdataset
```

Evaluate both Replogle test settings:

```bash
SPLITS=hf://datasets/arcinstitute/PIE_splits@396ab9563175ee887750c9eed7ccaea6f5fdbf50
for setting in test_seen test_unseen; do
  pie eval experiment_name=replogle_xdataset \
    run_dir=models/replogle_xdataset ckpt=best_auprc \
    split_path=$SPLITS/replogle_xdataset/$setting.json row_set=$setting
done
```

To write predictions for the unseen-perturbation test pairs:

```bash
SPLITS=hf://datasets/arcinstitute/PIE_splits@396ab9563175ee887750c9eed7ccaea6f5fdbf50
pie infer experiment_name=replogle_xdataset \
  run_dir=models/replogle_xdataset ckpt=best_auprc \
  rows_kind=split rows_path=$SPLITS/replogle_xdataset/test_unseen.json \
  output_path=models/replogle_xdataset/infer/test_unseen.parquet
```

Prediction parquet files contain one row per pair, arrays of the three outputs, and the gene
axis in file metadata. For a new context, prepare a `controls_only=true` dataset, provide a
query JSON file, and use the `preprocessed_dirs` override as described in the PIE documentation.
The original training assets remain necessary for response evidence. Runtime split files must
match the training split hash embedded in the checkpoint.

## License

The PIE model checkpoints and accompanying files are released under the
[Arc Research Institute PIE Model Non-Commercial License](MODEL_LICENSE.md) and are subject to
the [PIE Model Acceptable Use Policy](MODEL_ACCEPTABLE_USE_POLICY.md). The PIE code is licensed
separately under CC BY-NC-SA 4.0; see the [code repository](https://github.com/ArcInstitute/pie).
