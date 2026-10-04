# Per-cell similarity to external control references

This analysis compares **every one of the 18,400 non-targeting cells in each of A, B and C** with control-expression profiles from K562, HCT116, HEK293T, H1 and CD4. It produces 55,200 rows with five primary correlations per cell. This is an analysis of basal expression, not knockdown-response prediction or cell-line identification. Cells in different contexts are not paired: row 1 in A is not a biological match for row 1 in B.

## 1. Inspect existing inputs and source definitions

The scripts remain in the project root:

- `prepare_reference_controls.py`: load local controls and extract bounded remote control subsets.
- `analyze_cell_reference_correlations.py`: harmonize gene axes, calculate correlations for all competition cells, summarize and plot.
- `validate_cell_reference_correlations.py`: independently verify complete output axes, summaries, checksums and selected actual-cell correlations.
- `test_cell_reference_correlations.py`: check the numerical implementation against NumPy/SciPy and test deterministic selection and invalid inputs.

Existing `build_submission.py`, `pilot_orion.py` and `extract_cd4_reference.py` were read before implementing this analysis. The existing CD4 DE cache contains log fold changes, not unperturbed expression; it is not used as a basal profile. Instead we extract raw non-targeting counts from the authors' pseudobulk release.

Sources:

- K562: [Replogle et al. original processed release](https://api.figshare.com/v2/articles/20029387), `K562_gwps_raw_bulk_01.h5ad`. This release has a restricted measured gene axis; absent genes are never assumed to be zero.
- HCT116 and HEK293T: [Xaira X-Atlas-Orion](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion), pinned revision `53a5bc98d49247bcf967500292575c3d3602de31`, cached `Batch1` files and gene-token metadata. These are batch-limited profiles, not complete-atlas profiles.
- H1: [Arc's public 2025 train release](https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/train/adata_Training.h5ad?generation=1765904883947296). This older competition is an external reference for 2026; no 2026 perturbed measurements are used.
- CD4: [authors' source repository](https://github.com/emdann/GWT_perturbseq_analysis_2025), [data definitions](https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/data_sharing_readme.md), `GWCD4i.pseudobulk_merged.h5ad`. Its X matrix is summed raw UMI counts, not the DE-statistics matrix.

## 2. Select only non-targeting controls

K562 uses non-targeting core-control pseudobulk groups containing at least 30 cells. Each raw mean-expression row is multiplied by its cell membership before pooling; duplicate symbols are summed.

HCT116/HEK293T use all non-targeting cells in their cached first batch that pass the guide filter and have positive total counts. Extracted gene-expression sums must exactly match their annotated total counts.

H1 uses the first 32 non-targeting controls, in source row order, from each of the 48 batches in the 2025 train release: 1,536 controls, drawn from 38,176 available controls. It is a deterministic stratified sample, not a full release or a random sample.

CD4 uses the first 32 eligible non-targeting pseudobulks per donor/culture-condition/run stratum. Eligibility requires NTC annotation, positive cell membership, and passing minimum-cell and total-count checks. Its 384 selected pseudobulks cover 12 observed strata, represent 34,620 cells, and are drawn from 10,303 eligible pseudobulks. Raw row totals are checked against metadata. The primary CD4 profile pools all three conditions. Rest, Stim8hr and Stim48hr are also calculated separately as sensitivity comparisons.

The 32-per-stratum limit applies **only to external reference controls**. All 55,200 competition control cells are included without subsampling or expression-based exclusion.

HTTP requests use byte ranges and `If-Match` ETags. Each response must have the expected status, Content-Range and byte count. Cached range blocks have SHA-256 sidecars. Remote files are not fully downloaded. The default extraction budget is 2,048 MiB per source, including metadata ranges used. The HDF5 chunk-prefetch step only fetches storage chunks overlapping selected rows and is bounded by the same budget.

## 3. Harmonize genes and expression scales

Gene symbols are matched exactly, retaining the competition gene order. Duplicate reference symbols are summed. The primary analysis uses the **intersection of measured symbols across all five references and the competition**. No unmeasured gene is filled with zero. A `gene_axis_coverage.csv` records which genes each reference measures and which enter the common axis.

For each competition cell and each reference control profile:

1. Subset to this same common gene axis.
2. Divide each gene count by the total count on that axis and multiply by 1,000,000 (CPM).
3. Apply the natural-log transformation `log1p(CPM)`.
4. Calculate Pearson correlation across the common genes.

The reference is a pooled count profile. It is not the average of per-cell logged expression. K562 uses reconstructed count sums from means and membership. This weighting convention is consistent with using pooled counts for the other references.

All measured zeros are retained. Raw library size and detected-gene count in the output refer to the full original 18,533-gene competition matrix, not only the intersection.

## 4. Calculate sensitivity analyses

- Spearman correlation uses average ranks for ties, including tied zeros.
- Pair-specific Pearson correlations use the gene intersection for each reference separately. These columns have different gene axes and should not be used to rank references against each other.
- Positive-reference-gene Pearson correlations restrict the common axis to genes with positive pooled counts in all five references and renormalize that subset.
- Additional CD4 condition-specific profiles are compared on the primary common axis.
- Pooled A/B/C profiles are compared with each reference and with each other. Pooled correlations differ from the average or median of cell-level correlations.
- Per-guide mean correlations and associations with log library size help assess guide/depth dependence.

## 5. Outputs and interpretation

Outputs are written to `analysis/cell_reference_correlations/`:

| File | Contents |
|---|---|
| `per_cell_correlations.csv` | All 55,200 competition cells; original cell ID and guide ID; full-axis library/detection metadata; Pearson, Spearman and sensitivity comparisons; closest reference and the margin to second place |
| `context_reference_summary.csv` | Mean, median, spread, 5th/95th percentiles, closest-reference fractions and depth associations for each context/reference |
| `cd4_condition_summary.csv` | Condition-specific CD4 median and 5th/95th percentile correlations |
| `pooled_context_correlations.csv` | Correlations between pooled competition controls and references |
| `context_profile_correlations.csv` | A/B/C pooled-profile correlations on the same common gene axis |
| `reference_profile_correlations.csv` | Pairwise correlations among the five pooled reference profiles |
| `gene_axis_coverage.csv` | Full competition axis and reference coverage flags |
| `reference_control_profiles_common_genes.csv` | The actual reference count profiles used on the common axis |
| `reference_selected_controls_*.csv` | Exact external source rows selected, plus their available metadata |
| `analysis_report.json` | Source identity, input/output checksums, dimensions, versions, validation and limitations |
| `correlation_overview.png` | Per-context median/spread and closest-reference frequencies |

A closest reference means **closest among these five profiles under this normalization**. High correlation can reflect common housekeeping genes, shared zeros, cell state, count depth or platform effects. It does not establish the identity of A/B/C, and it does not prove that the closest reference has the most transferable CRISPRi responses. Reference batch and condition selection are additional limitations.

## 6. Reproduce in Windows Command Prompt

Place the original `context_A.h5ad`, `context_B.h5ad` and `context_C.h5ad` validation controls in the project root alongside the tracked manifest and gene/target lists. The large input matrices are not included in Git. On a fresh Windows machine, create the environment first with `py -3.14 -m venv .competition_env`.

Run from the project root using the pinned environment:

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m pip install -r requirements-competition.txt
.\.competition_env\Scripts\python.exe download_reference.py
.\.competition_env\Scripts\python.exe pilot_orion.py
.\.competition_env\Scripts\python.exe prepare_reference_controls.py
.\.competition_env\Scripts\python.exe -m unittest test_cell_reference_correlations -v
.\.competition_env\Scripts\python.exe analyze_cell_reference_correlations.py
.\.competition_env\Scripts\python.exe validate_cell_reference_correlations.py
```

`pilot_orion.py` also reproduces the earlier audit outputs; it is needed here only if the Orion cache is absent. `download_reference.py` and `pilot_orion.py` reuse verified local downloads. The default reference extraction uses 32 controls per stratum; changing this setting changes the reference and should be treated as a new sensitivity experiment.

To inspect remote schemas without extracting counts:

```cmd
.\.competition_env\Scripts\python.exe prepare_reference_controls.py --inspect H1
.\.competition_env\Scripts\python.exe prepare_reference_controls.py --inspect CD4
```

To repeat just one extraction:

```cmd
.\.competition_env\Scripts\python.exe prepare_reference_controls.py --source H1 --controls-per-stratum 32 --budget-mib 2048
.\.competition_env\Scripts\python.exe prepare_reference_controls.py --source CD4 --controls-per-stratum 32 --budget-mib 2048
```

The range cache and prepared source profiles live under `reference/control_profiles/` and are excluded from Git. Derived numeric results and provenance are retained in the analysis output directory. The analysis creates no predictions and performs no leaderboard upload.


## Relationship to the earlier pilot

`analysis/model_audit/reference_context_correlations_common_axis.csv` is an earlier three-reference, context-average pilot. It is retained as historical evidence. The new analysis adds H1/CD4, changes to a five-reference common axis, uses pooled control counts for references and pooled-context summaries, and evaluates every individual competition cell. Its numbers are therefore not direct repeats of that pilot; a change in correlation can reflect gene-axis or aggregation choices.

The focused tests cover independent Pearson/Spearman oracles (including tied zeros and constant vectors), CPM scale invariance, invalid counts, deterministic stratum selection, range byte budgets and corrupted cache detection. The full-data validator checks all cell identifiers and guides, summary medians/winner counts, every recorded output CSV hash and 45 cell/reference pairs calculated directly from original counts.

## Completed result

The completed run used 7,447 common genes. Median cell-to-reference Pearson correlations are highest for CD4 in A (0.4821), and H1 in B (0.4635) and C (0.4766). These references win for 18,399/18,400 A cells, 18,393/18,400 B cells and 18,397/18,400 C cells. Spearman and the positive-reference-gene sensitivity retain the same leading reference at the context-median level. Pooled-profile results also retain these leaders. The pooled B/C correlation is 0.7210, A/B is 0.6293, and A/C is 0.5073 on the same gene axis.

All 55,200 output cell identifiers and guides, CSV checksums and summary medians/winner counts passed independent validation. Forty-five actual-cell/reference pairs matched independent NumPy/SciPy calculations. The validator explicitly converts original float32 storage to float64 before normalization, avoiding small rounding differences. Three focused numerical/extraction tests pass.

Full tables are regenerated in [CELL_REFERENCE_CORRELATION_RESULTS.md](CELL_REFERENCE_CORRELATION_RESULTS.md). These are control-expression associations and should not be used as cell identity assignments or automatic perturbation-transfer weights.
