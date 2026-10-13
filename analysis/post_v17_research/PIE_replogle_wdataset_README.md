---
library_name: pytorch
license: other
license_name: pie-model-license
license_link: https://huggingface.co/arcinstitute/PIE_replogle_wdataset/blob/main/MODEL_LICENSE.md
tags:
  - biology
  - gene-expression
  - perturbation-prediction
  - perceiver
  - cross-validation
datasets:
  - arcinstitute/PIE_replogle_nadig_essential
  - arcinstitute/PIE_sources
---

# PIE: within-dataset generalization on Replogle

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

This repository groups four independently trained models from the `replogle_wdataset`
experiment. Each model holds out one Replogle cell line and trains and validates on the other
three. The four models are separate folds; the reported results do not use an ensemble.

## Checkpoints and files

**The folder name identifies the held-out cell line**, rather than a training cell line.
For example, `k562/best_auprc.ckpt` was trained and validated on HepG2, Jurkat, and RPE1.

| Directory | Held-out cell line | Training and validation cell lines |
|---|---|---|
| `hepg2/` | HepG2 | Jurkat, K562, RPE1 |
| `jurkat/` | Jurkat | HepG2, K562, RPE1 |
| `k562/` | K562 | HepG2, Jurkat, RPE1 |
| `rpe1/` | RPE1 | HepG2, Jurkat, K562 |

```text
README.md
hepg2/
  best_auprc.ckpt
  config.yaml
  data_stats.json
  eval/
    unseen_ctx/{metrics_best_auprc.csv,granular_best_auprc.csv}
    unseen_pert/{metrics_best_auprc.csv,granular_best_auprc.csv}
    unseen_ctx_pert/{metrics_best_auprc.csv,granular_best_auprc.csv}
jurkat/  # same layout
k562/    # same layout
rpe1/    # same layout
```

`best_auprc.ckpt` is the Lightning checkpoint with the highest validation `binary_auprc`.
The checkpoint embeds its model configuration, pinned asset references, and fitted data
statistics; the adjacent YAML and JSON files make those settings inspectable.

## Data and training

Training uses the preprocessed
[Replogle-Nadig essential-gene dataset](https://huggingface.co/datasets/arcinstitute/PIE_replogle_nadig_essential)
at revision `20c9faef76fc96fdc809871340bd93499413dfe7`. Its measured gene axis has 6,642 genes.
Canonical train, validation, and test splits live in
[PIE splits](https://huggingface.co/datasets/arcinstitute/PIE_splits) under `replogle_wdataset/`.
Each fold uses the same train and validation splits for all three test
settings.

Knowledge sources are `esm2`, `ncbi_text`, `string_space`, `depmap_gene_effect`, `context_text`,
and `perturbation_text`, plus `gene_text` for gene queries. All come from
[PIE sources](https://huggingface.co/datasets/arcinstitute/PIE_sources) at revision
`cb1aaa4e7655605bdc70a9bd77bbd62016b8c7d7`.

## Evaluation

All results use `best_auprc.ckpt` and the canonical test splits.

| Setting | Test cell contexts | Test perturbations |
|---|---|---|
| `unseen_ctx` | Held-out cell line | Seen in training |
| `unseen_pert` | The three training cell lines | Unseen in training |
| `unseen_ctx_pert` | Held-out cell line | Unseen in training |

Values below are **mean ± one sample standard deviation across four folds** (`ddof=1`).

| Setting | AUPRC | Jaccard | Direction match | Spearman DE count | Spearman LFC | L1 discrimination |
|---|---|---|---|---|---|---|
| `unseen_ctx` | 0.4155 ± 0.0286 | 0.1771 ± 0.0128 | 0.8235 ± 0.0206 | 0.7389 ± 0.0395 | 0.6568 ± 0.0411 | 0.7701 ± 0.0604 |
| `unseen_pert` | 0.3048 ± 0.0256 | 0.1670 ± 0.0114 | 0.7775 ± 0.0135 | 0.5420 ± 0.0350 | 0.5505 ± 0.0201 | 0.7142 ± 0.0106 |
| `unseen_ctx_pert` | 0.2366 ± 0.0402 | 0.1160 ± 0.0198 | 0.7257 ± 0.0197 | 0.4658 ± 0.0717 | 0.4475 ± 0.0469 | 0.6526 ± 0.0429 |

Metric names in the CSVs are `binary_auprc`, `sig_jaccard`, `direction_match`, `spearman_nsig`,
`spearman_lfc`, and `discrimination_score_l1`. Higher is better for all six.

## Download and use

Install PIE with `pip install arc-pie` (or, from a checkout of the PIE code repository,
`uv sync --frozen` for the pinned environment). Set the writable data, runs and cache
directories in `common.sh` (copy `common.sh.example` from the code repository):

```bash
pip install arc-pie
hf auth login

# Download the K562-held-out model and its inspectable configuration.
hf download arcinstitute/PIE_replogle_wdataset \
  k562/best_auprc.ckpt k562/config.yaml k562/data_stats.json \
  --local-dir models/replogle_wdataset
```

Evaluate the K562-held-out model on each generalization setting:

```bash
SPLITS=hf://datasets/arcinstitute/PIE_splits@396ab9563175ee887750c9eed7ccaea6f5fdbf50
for setting in unseen_ctx unseen_pert unseen_ctx_pert; do
  pie eval experiment_name=replogle_wdataset/k562 \
    run_dir=models/replogle_wdataset/k562 ckpt=best_auprc \
    split_path=$SPLITS/replogle_wdataset/$setting/k562/test.json \
    row_set=$setting
done
```

To write predictions for the held-out-context test pairs:

```bash
SPLITS=hf://datasets/arcinstitute/PIE_splits@396ab9563175ee887750c9eed7ccaea6f5fdbf50
pie infer experiment_name=replogle_wdataset/k562 \
  run_dir=models/replogle_wdataset/k562 ckpt=best_auprc \
  rows_kind=split rows_path=$SPLITS/replogle_wdataset/unseen_ctx/k562/test.json \
  output_path=models/replogle_wdataset/k562/infer/unseen_ctx.parquet
```

Replace `k562` with another fold to use its checkpoint and splits. Prediction parquet files
contain one row per pair, arrays of the three outputs, and the gene axis in file metadata.
For a new context, prepare a `controls_only=true` dataset, provide a query JSON file, and use
the `preprocessed_dirs` override as described in the PIE documentation. The original training
assets remain necessary for response evidence. Runtime split files must match the training
split hash embedded in the checkpoint.

## License

The PIE model checkpoints and accompanying files are released under the
[Arc Research Institute PIE Model Non-Commercial License](MODEL_LICENSE.md) and are subject to
the [PIE Model Acceptable Use Policy](MODEL_ACCEPTABLE_USE_POLICY.md). The PIE code is licensed
separately under CC BY-NC-SA 4.0; see the [code repository](https://github.com/ArcInstitute/pie).
