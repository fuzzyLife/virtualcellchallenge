# V2 submission: donor-calibrated CD4 hybrid

Generated for the validation panel `vcc2026-val-1` on 2026-10-04. The candidate is separate from the first submission; it is published with overall score -0.22332546152913602 (entry `ENmII1VynYrgpB2KICPC`).

## Model and evidence

The original audit found severe underprediction and weak separation of K562 effects from source-control noise. V2 uses a statistically modeled public CD4 Perturb-seq reference for eligible targets, with the existing K562 response model as fallback. It uses source log2 fold changes, not source expression counts, and generates raw integer counts from the challenge's own control cells.

- CD4 contributes to 223 targets; 61 use K562 fallback; 16 have only a modeled on-target knockdown.
- Combined reference coverage is 284/300 targets and 11,163/18,533 genes. Coverage describes measured source axes, not the number of significant predicted genes.
- CD4 rows require at least 100 cells, significant on-target knockdown, two guides, and no reported adjacent-gene, distal off-target or low-target-expression flags.
- Settings were selected using two disjoint donor pairs: `CE0006864_CE0008162` and `CE0008678_CE0010866`. There were 138 eligible shared target-condition contrasts covering 72 targets and 10,273 genes.
- A predeclared grid compared source FDR cutoffs 0.01/0.05/0.1/0.2 and log2-effect shrinkage 0/0.3/0.5/0.7/1.0. Bidirectional donor prediction error excludes the target gene and clips log2 fold changes to +/-3. The selected cutoff is 0.2 and shrinkage is 0.3, with relative squared error 0.986857 against the zero-response anchor (about 1.3% lower).
- Eligible CD4 conditions are averaged with square-root cell-count weights. No cell-line identity or condition assignment is inferred for A/B/C.
- K562 supplies genes unmeasured in CD4, or targets without an eligible CD4 condition. It retains the original z=1.96 source-noise threshold with a fixed 0.7 transfer coefficient. This fallback coefficient is exploratory and is not calibrated by the CD4 donor test.
- Upward count simulation can activate genes with zero counts in sampled cells, using context-wide gene means and cell library factors. Downward responses use binomial thinning; target genes retain the original assumed 90% knockdown.

The `--transfer-shrinkage 0.7` CLI setting applies to K562 fallback in this model. CD4 shrinkage is separately selected from donors and recorded in `reference_model.json`.

There are 18,977 modeled downstream gene-target effects above absolute log2FC 0.1, compared with 1,871 in the original model. This increase does not establish accuracy. The donor error is an exploratory source-domain result, enriched for larger responses, with no third validation split; it is neither a cross-cell-line result nor a competition score. Unsupported genes receive no downstream effect, and all three contexts still share the source effect multipliers.

## Public data extraction and provenance

Source: [Zhu/Dann et al. primary CD4 Perturb-seq resource](https://doi.org/10.1016/j.cell.2026.08.002), with the [authors' public data description](https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/data_sharing_readme.md). Public external data is permitted according to the [2026 competition announcement](https://arcinstitute.org/news/virtual-cell-challenge-2026).

`extract_cd4_reference.py` downloads only relevant contiguous HDF5 byte ranges. Each request requires the recorded object size, ETag and exact Content-Range. The extracted merged subset has 854 rows; each selected donor subset has 190 rows. The three DE layers total 325,885,488 downloaded bytes, plus bounded metadata. Full 16.8 GB source files and the 44.6 GB pseudobulk release are not downloaded.

Checksums and source object identities are in `analysis/model_audit/cd4_extraction.json`. Compressed subsets and metadata CSVs remain locally under `reference/cd4_probe/`, excluded from Git. Source definitions distinguish log2 fold change (`log_fc`), adjusted p-value (`adj_p_value`) and standard error (`lfcSE`); on-target summary statistics are not substituted for log fold changes.

## Reproduce (Windows Command Prompt)

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe extract_cd4_reference.py
.\.competition_env\Scripts\python.exe -m unittest test_submission test_cd4_model -v
.\.competition_env\Scripts\python.exe build_submission.py --output-dir experiments/v2_cd4_hybrid --reference-kind cd4-hybrid --transfer-shrinkage 0.7 --noise-z 1.96 --upward-mode zero-aware
.\.competition_env\Scripts\python.exe validate_v2.py
```

An existing completed output requires intentional `--force` to rebuild. Original defaults and the first root-level prediction are retained.

## Submit (Windows Command Prompt)

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v2_cd4_hybrid\prediction.vcc -m "Donor-calibrated CD4 hybrid v2" --wait
```

The existing keychain login can be reused. Upload and scoring are separate from generation; the generator does not submit the package. Results should be recorded separately from the first entry `QQnFEBGlOk2pFFwkgtJL`.

## Completed validation and response check

All 360,000 cells and 900 groups pass the full raw-count validator: exactly 400 cells per target/context, exact gene order and labels, positive stored integer counts, 2,083,968,482 nonzero entries, and maximum cell library size 52,700 (limit 1,000,000). All 15 distinct pipeline/source-model tests pass. The original package SHA-256 was freshly verified unchanged. Downloaded source subset checksums also match their extraction records.

The paired diagnostic checks 12 deliberately selected targets spanning source coverage and effect strength, in all contexts. SMARCA5 produces a much larger response; weaker targets can remain silent. This is not a random panel estimate and does not establish true DE recall or direction accuracy.

Significant downstream genes against the provided controls (original -> v2):

| Target | A | B | C |
|---|---:|---:|---:|
| TMEM41B | 0 -> 3 | 0 -> 4 | 0 -> 1 |
| R3HCC1L | 0 -> 9 | 0 -> 6 | 0 -> 2 |
| SMARCA5 | 4 -> 558 | 2 -> 549 | 3 -> 506 |
| ABCD1 | 0 -> 0 | 0 -> 0 | 0 -> 0 |

Detailed records: `experiments/v2_cd4_hybrid/analysis/paired_response_check.csv`, `prediction_validation.json`, `donor_validation.json`, and `provenance_check.json`. Package readiness is confirmed only by a completed official prep run and the generated `submission_ready.json`.

## Upload-ready artifact

Official prep and package validation passed. File: `experiments/v2_cd4_hybrid/prediction.vcc` (3,307,438,080 bytes; 3.307 GB).

SHA-256: `faa73819d4638e4179d9d2b8dac726a1516cfce5d85a867a87c1e152c1b39656`.

The package was uploaded successfully and the entry is published. Its checksum differs from the original package. Full readiness and official package metadata are in `experiments/v2_cd4_hybrid/analysis/submission_ready.json`.

## Published leaderboard result

Entry `ENmII1VynYrgpB2KICPC` was published on 2026-10-04. Rank was 1105 in the upload output and 1104 at the subsequent read-only API check; rank changes as other entries change. The panel and anchor set match the first submission.

| Scaled component | First | V2 | Change |
|---|---:|---:|---:|
| overall | -0.276290 | -0.223325 | +0.052964 |
| pds | 0.019644 | 0.103042 | +0.083398 |
| mse | 0.000000 | 0.000000 | +0.000000 |
| nmae | 0.016006 | 0.032352 | +0.016347 |
| fid | -1.665331 | -1.470413 | +0.194919 |
| reach | 0.033874 | 0.041380 | +0.007505 |
| jac | -0.061932 | -0.046314 | +0.015618 |

The overall gain is +0.052964. All five non-expression scaled components improve. Raw expression error also improves, from 1.023018 to 1.014294, but remains on the zero-score floor. Raw discrimination rises from 0.509915 to 0.547923, still fairly close to chance (0.5). Raw direction-fidelity yield rises from 0.017572 to 0.074675; this measures correct significant calls relative to the larger response budget, not the fraction of predicted signs that are correct. Raw significant-set Jaccard rises from 0.007733 to 0.013545.

V2 is a measurable improvement on this validation panel and remains below the zero-score anchor. The comparison changes reference data, filtering, coverage, transfer and simulation together, so it cannot attribute the gain to one component. Direction fidelity remains the dominant negative term. Next improvements should address context-dependent response transfer and correctly directed response coverage, while checking source specificity. Aggregate metrics cannot separate missed responses from wrong predicted directions.

Records: `leaderboard_result.json`, `published_api_status.json` and `leaderboard_comparison.csv` in the experiment analysis directory. The first submission's historical result is retained.
