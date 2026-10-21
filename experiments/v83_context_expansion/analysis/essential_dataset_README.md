---
pretty_name: PIE Replogle-Nadig essential Perturb-seq (4 cell lines)
tags:
  - biology
  - single-cell
  - perturb-seq
  - crispri
  - differential-expression
size_categories:
  - 10M<n<100M
configs:
  - config_name: k562
    data_files: de/K562.parquet
  - config_name: rpe1
    data_files: de/RPE1.parquet
  - config_name: hepg2
    data_files: de/hepg2.parquet
  - config_name: jurkat
    data_files: de/jurkat.parquet
---

# PIE: Replogle-Nadig essential Perturb-seq

Differential-expression (DE) tables and the preprocessed dir for
[PIE](https://github.com/ArcInstitute/pie), built from the CRISPRi essential-gene Perturb-seq screens
in K562 and RPE1 (Replogle et al., 2022) and HepG2 and Jurkat (Nadig et al., 2025). This dataset is
the test set of the PIE `replogle_xdataset` experiment and the train and test data of
`replogle_wdataset`.

Datasets and knowledge sources: [PIE collection](https://huggingface.co/collections/arcinstitute/pie-6abc354d475081a4a7cb4492)

## Contents

```
de/<cell_line>.parquet   per-cell-line DE, one row per (perturbation, gene)
preprocessed/            PIE preprocessed dir, with the context map contexts.yaml
```

| Cell line | Cells   | Perturbed genes | DE rows   |
|-----------|---------|-----------------|-----------|
| K562      | 179,753 | 1,314           | 8,727,588 |
| RPE1      | 167,372 | 1,453           | 9,650,826 |
| HepG2     | 90,826  | 1,253           | 8,322,426 |
| Jurkat    | 173,260 | 1,439           | 9,557,838 |

Controls are the `non-targeting` cells of each cell line.

### `de/`

One parquet per cell line, from the `pie process` `de` stage (gpudge 0.7.0 on raw counts, grouped
by `target_gene` against `non-targeting`, CP10k normalization, genes kept at ≥ 5 CPM).

| Column             | Description                                        |
|--------------------|----------------------------------------------------|
| `target`           | perturbed gene                                     |
| `feature`          | measured gene                                      |
| `target_mean`      | mean expression in perturbed cells                 |
| `ref_mean`         | mean expression in control cells                   |
| `target_ncells`    | number of perturbed cells                          |
| `ref_ncells`       | number of control cells                            |
| `log2_fold_change` | log2 fold change, perturbed vs. control            |
| `p_value`          | Mann-Whitney U test p-value                        |
| `Ueffect`          | U-statistic effect size                            |
| `p_adj`            | Benjamini-Hochberg adjusted p-value                |

### `preprocessed/`

The PIE preprocessed dir (`pie prep` output, format version 1) that `pie sources`, `pie train`,
`pie eval` and `pie infer` read: N = 5,459 rows, one per (cell line, perturbation) pair, over
G = 6,642 genes and C = 4 cell lines.

| File               | Shape, dtype      | Description                                                    |
|--------------------|-------------------|----------------------------------------------------------------|
| `meta.json`        |                   | gene axis, cell line and perturbation vocabularies, checksums, perturbation Ensembl ids |
| `contexts.yaml`    |                   | context -> Cellosaurus accession map; `pie sources` reads it for `context_text` |
| `fold_changes.npy` | (N, G) float32    | linear fold change, 0 where untested                           |
| `fdr.npy`          | (N, G) float32    | adjusted p-value, 1 where untested                             |
| `tested.npy`       | (N, G) bool       | whether the gene was tested for the pair                       |
| `lfc_true.npy`     | (N, G) float64    | log2 fold change, NaN where untested                           |
| `delta_p.npy`      | (N, G) float32    | pseudobulk log1p expression of perturbed cells minus control mean |
| `ctrl_means.npy`   | (C, G) float32    | mean log1p expression of control cells per cell line           |
| `ctx_ids.npy`      | (N,) int32        | cell line of each row, index into `context_to_id`              |
| `pert_ids.npy`     | (N,) int32        | perturbation of each row, index into `pert_to_id`              |

## Usage

The `replogle_wdataset` and `replogle_xdataset` experiment configs pin this repo by commit, so `pie
train` and `pie eval` download `preprocessed/` automatically. In your own config, use
`hf://datasets/arcinstitute/PIE_replogle_nadig_essential@<commit>/preprocessed` in
`data.preprocessed_dirs`, or download a copy:

```bash
hf download arcinstitute/PIE_replogle_nadig_essential --repo-type dataset \
  --local-dir "$PIE_DATA_ROOT/replogle"
```
