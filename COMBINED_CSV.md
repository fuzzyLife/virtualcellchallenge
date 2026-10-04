# Combined competition control CSV

`competition_controls_combined.csv` contains all supplied competition control expression data in one file. It combines `context_A.h5ad`, `context_B.h5ad`, `context_C.h5ad`, `gene_names.csv`, `pert_counts.csv` and `manifest.json`. External reference experiments and generated predictions are excluded.

## Layout

There is one header row, 18,533 gene rows, 55,200 numeric control-cell columns and 31 annotation columns (55,231 columns total). Gene order matches the original `gene_names.csv`; cells retain their original order within A, then B, then C. Counts are unchanged raw integers, including all original zeros. No cells or genes are sampled, filtered or normalized.

Example cell column:

```text
context_A_col00001__cell_id=A_ctrl_000000__ntc_id=non-targeting-11__total_counts=<integer>__detected_genes=<integer>
```

The prefix identifies context and one-based cell-column position. The suffix records the original cell ID, original control-guide ID, total raw counts and detected-gene count. Identifier values are URL-encoded where needed and can be recovered with `urllib.parse.unquote`. Metadata are kept in the header so every cell-count column remains numeric and the CSV has no mixed metadata rows.

## Annotation columns

- `gene_name`, `gene_index_1based`: original gene symbol and feature position.
- `is_requested_knockdown_target`: 1 for the 300 requested targets, otherwise 0. This is a gene-list annotation, not an observed perturbed-cell label.
- `requested_target_index_1based`: original row position in `pert_counts.csv`; 0 for other measured genes.
- `requested_cells_per_target_context`: 400 for requested targets, otherwise 0. This describes required predictions, not additional measured cells.
- `season`, `partition`, `panel_id`, `source_url`, `count_unit`: panel provenance and measurement unit.
- `observed_target_gene`: always `non-targeting` in the supplied controls.
- `has_observed_knockdown_responses`: always 0.
- `panel_gene_count`, `panel_target_count`, `panel_required_cells_per_target`, `panel_required_prediction_cells_total`: original panel dimensions and requested prediction counts.
- For each A/B/C: `context_*_control_cells`, `context_*_control_guides`, `context_*_ground_truth_cells_declared`, `context_*_gene_total_raw_count`, `context_*_gene_detected_control_cells`.

Declared ground-truth cell counts are manifest metadata only. Hidden response measurements are not present in this export.

## Training use

This is a complete **control-expression dataset**, suitable for learning basal expression patterns or control-cell representations. It contains no measured knockdown outcomes, so it cannot by itself train a supervised perturbation-response model. Target flags must not be treated as per-cell intervention labels. The usual cell-by-gene machine-learning matrix is the transpose of the context-prefixed count columns; annotation columns should be kept separately.

Use Python/R for the full file. Its 55,231 columns exceed Excel's column capacity. Stream gene rows in chunks or select columns rather than loading the complete mixed table as strings.

```python
import pandas as pd

path = 'competition_controls_combined.csv'
header = pd.read_csv(path, nrows=0).columns
count_columns = [c for c in header if c.startswith(('context_A_col', 'context_B_col', 'context_C_col'))]

# Read raw counts as int32 in gene blocks; all other columns are annotations.
for block in pd.read_csv(path, chunksize=64, dtype={c: 'int32' for c in count_columns}):
    gene_names = block['gene_name'].to_numpy()
    counts_for_gene_block = block[count_columns].to_numpy()  # genes x cells
    # Process this block without retaining the full dataset in memory.
```

## Reproduce and verify

```cmd
.\.competition_env\Scripts\python.exe -m unittest test_export_competition_csv -v
.\.competition_env\Scripts\python.exe export_competition_csv.py
```

An existing CSV requires `--force` for intentional replacement. `--phase extract`, `--phase write` and `--phase verify` allow separate extraction, generation and verification. The disposable cache is under `analysis/competition_csv_cache/`.

Verification reads the complete CSV back, checks all headers, annotation fields and gene order, verifies a SHA-256 over every count in canonical gene-major int32 order, and independently checks every cell total and every context/gene total. The test fixture compares against independent dense matrix transposes and confirms that a corrupted count is detected. Results and input-file checksums are recorded in `analysis/combined_csv_export.json`.

The CSV and transpose cache are excluded from Git because of their size. The exporter, test and documentation remain reproducible project code.

## Completed export

The file contains 18,533 gene rows and 55,231 columns, including all 55,200 control cells and 300 target flags. Size: 2,076,031,902 bytes (2.076 GB). The complete numeric and metadata read-back verification passed.

CSV SHA-256: `595a7f693b4638627059abbe82c58b468d39f95b4733c09430682c12f8892013`.
