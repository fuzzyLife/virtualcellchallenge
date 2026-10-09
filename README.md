# Virtual Cell Challenge 2026

This project analyses the validation bundle for the [Virtual Cell Challenge](https://virtualcellchallenge.org/), generates raw single-cell CRISPRi response predictions, and packages them using the official VCC CLI. Nine submissions (first, V2, V4, V5, V6, V7, V9, V10 and V11) have been published. V11 is the best recorded validation submission, with overall score **0.137003974**, improving on V10 by **0.020249853**. Expression accuracy and biological dispersion remain weaknesses. V3 remains an unsubmitted candidate. See [V11_LEADERBOARD_RESULTS.md](V11_LEADERBOARD_RESULTS.md) for the current result and [V11_SUBMISSION.md](V11_SUBMISSION.md) for reproducibility.

The latest and best published submission is V11 (`bmIAGveuEvFiYQTvz0Ma`), overall **0.137003974**. See [V11_LEADERBOARD_RESULTS.md](V11_LEADERBOARD_RESULTS.md).

Reviewed the paper, code and cached embeddings in `cross_modal_foundation_model`. The root script `review_cross_modal_model.py` records provenance, exact gene coverage and synthetic checks. See [CROSS_MODAL_MODEL_REVIEW.md](CROSS_MODAL_MODEL_REVIEW.md): RNA-only Tahoe embeddings cover 235/300 targets, while the advertised proteomics checkpoint is currently absent. A context-aware response-head pilot is recommended; no submission was changed.

Implemented and evaluated four interventional response heads on 889 public contrasts (284 targets), with whole-lineage and target holdouts. Tested 180 matched generated-cell comparisons on nine larger H1 target groups. The embedding/PCA blends worsened DE fidelity, overlap and corrected expression error, so V4 was retained and no new upload package was generated. See [EMBEDDING_RESPONSE_PILOT.md](EMBEDDING_RESPONSE_PILOT.md) for scripts, results and reproduction.

V7 is published as `k0pkYuJRcw3C0U8K0NgL`, overall **-0.051386751** versus V4 **-0.050673910**. It uses exactly V4 response effects with new count/control sampling. Expanded HCT116/HEK293T sources and neural/ridge alternatives were tested, but none passed the calibration guards. Official packaging and the independent full-payload audit passed. V4 was best at V7 publication. See [V7_LEADERBOARD_RESULTS.md](V7_LEADERBOARD_RESULTS.md), [V7_SUBMISSION.md](V7_SUBMISSION.md) and [V7_MODEL_STRATEGY.md](V7_MODEL_STRATEGY.md).

A fresh public-method review and a reproducible evaluation audit found that our historical local PDS used target-specific control origins. Identical predictions scored 0.854545-1.000000 locally instead of the required 0.5. See [VCC2026_RESEARCH_AND_GAPS.md](VCC2026_RESEARCH_AND_GAPS.md), [TRANSFER_VALIDATION_AUDIT.md](TRANSFER_VALIDATION_AUDIT.md) and root script `audit_transfer_validation.py`. This audit preserves prior submissions; the V8/V9 pipelines below implement corrected validation before selecting new count generators.

<!-- v8-status:begin -->

V8 is retained as a verified diagnostic pilot at `experiments/v8_corrected_transfer/prediction.vcc`, with official prep and independent payload verification passed. Candidate `baseline_dual` uses `v4_a100` and `dual`; public confirmation proxy improvement accepted: False. It has not been submitted. See [V8_SUBMISSION.md](V8_SUBMISSION.md) for the corrected validation, neutral controls and tradeoffs. V9 is published below; V11 is now the best published reference.

<!-- v8-status:end -->


<!-- V9_VERIFIED_SUBMISSION -->

## V9 published result

V9 is published as `QBpef7kQNiFa4IFSxNi5`, overall **-0.070447281**. Official prep and independent full-payload audit passed. Its 50% count calibration improved five scaled metrics versus V4 but worsened fidelity enough to lower overall performance. V4 biological effects are unchanged. Public confirmation acceptance did not predict leaderboard improvement; its subset fidelity guard was too permissive. V4 remained best at V9 publication; V11 is now best. See [V9_LEADERBOARD_RESULTS.md](V9_LEADERBOARD_RESULTS.md) and [V9_SUBMISSION.md](V9_SUBMISSION.md). V8 remains an unsubmitted diagnostic pilot.


<!-- post-v11-research:begin -->

## Post-V11 research and tests

Reviewed ten relevant scientific papers/methods and separately checked official data and participant reports. Collected the VIPerturb manifest (300/300 targets; 124 QC-filtered; median 47 cells), pinned Reactome annotations and RNA-fingerprinting code, and ran 48 control-distribution comparisons plus five-family transfer diagnostics. A gene-specific gamma-Poisson pilot reduces average null DE calls from 1,711 to 84 against 2,000 references, but retains context-specific bias and has no perturbation/leaderboard validation. Twelve unused count-confirmation targets are frozen for the next experiments. V11 remains the best submitted model. See [POST_V11_RESEARCH.md](POST_V11_RESEARCH.md) for papers, findings, limitations and reproduction.

<!-- post-v11-research:end -->

## Completed work

- Reviewed and hashed all 146 original project files, including the manifest, gene/target lists, three control matrices, saved competition pages and their assets.
- Extracted the saved HTML and all 61 PDF pages. The PDF concerns bacteriophage genome generation and was excluded from the model.
- Checked every stored control count, gene order, context label and guide count. Analysed per-gene expression, detection, library sizes, guide balance and differences between contexts.
- Rejected a public essential-gene reference because it had zero target overlap with this panel.
- Downloaded and checksum-verified the original Replogle et al. genome-wide K562 pseudobulk reference. Retained usable responses for 260/300 targets and 7,681/18,533 challenge genes.
- Implemented the root-level prediction pipeline, source downloader and five focused tests. Fixed a quantile-bin sampling overlap discovered during the initial full run.
- Generated and independently checked the complete prediction, passed official `vcc prep`, and created the upload package.
- Uploaded `prediction.vcc` as **Public CRISPR transfer**. The competition published entry **QQnFEBGlOk2pFFwkgtJL**.

The detailed review, modeling assumptions and limitations are in [COMPETITION_ANALYSIS.md](COMPETITION_ANALYSIS.md). No challenge perturbation ground truth was available locally.

## Competition requirements

This is the **validation** partition, panel `vcc2026-val-1`, with contexts **A/B/C**. Each context supplies 18,400 non-targeting controls from 46 guides.

| Requirement | Value |
|---|---|
| Targets | Exactly the 300 symbols in `pert_counts.csv`, in each context |
| Cells per target/context | 400 |
| Prediction shape | 360,000 cells x 18,533 genes |
| Gene order | Exact order in `gene_names.csv` (header: `gene_name`) |
| Cell annotations | `obs.target_gene` and `obs.context`, retaining A/B/C labels |
| Expression | Finite nonnegative whole raw counts in sparse CSR form |
| Controls in submission | None |
| Maximum counts per cell | 1,000,000 |
| Maximum stored entries | 4,750,000,000 |
| Upload file | One `prediction.vcc` covering all three contexts |

The target gene itself is excluded from all six scored metrics; the discrimination metric excludes all 300 panel target genes. Predicting only the direct knockdown therefore cannot solve the evaluated task.

The final round is separate: the supplied documentation schedules new contexts D/E/F and a new target panel for October 22, 2026, with submissions due November 5, 2026. This A/B/C artifact is a validation-round submission.

## Model

`build_submission.py` pools effective K562 source constructs (at least 30 cells and at least 50% source knockdown), aligns gene symbols and sums duplicate source symbols. It computes source log2 response ratios with a 1-CPM pseudocount, soft-thresholds sampling noise, clips effects, and shrinks estimates with limited source cells.

The pipeline transfers **30% of the estimated log2 effect** onto each context's own controls. This fixed shrinkage is an assumption, not a validated cross-context fit. The other 40 targets and genes outside the source axis retain control-derived downstream expression. A separate 90% on-target reduction is also an assumption; the target is excluded from scoring.

For each target/context, the pipeline samples 400 distinct control cells, balanced across guides and stratified by library size. It applies binomial thinning for decreases and Poisson increments for increases, then streams integer CSR blocks to H5AD with 64-bit global pointers. The default random seed is **20261003**.

Different promoter/transcript constructs were compared only as an exploratory source agreement diagnostic. The legacy filename `reference_cross_validation.csv` does **not** establish cross-context validation and was not used to tune the transfer setting.

## Validation and generated artifacts

| Artifact/check | Completed result |
|---|---|
| `prediction.h5ad` | 2,473,492,232 bytes |
| `prediction.vcc` | 3,306,905,600 bytes (3.31 GB) |
| Groups | 900; exactly 400 cells each |
| Context sizes | 120,000 cells each in A, B and C |
| Stored entries | 2,083,470,063 |
| Maximum predicted cell total | 52,457 counts |
| Gene order and labels | Exact matches |
| Stored values | Positive integers; explicit zeros removed |
| Official packaging | Passed with `vcc-cli==0.2.2` |
| Tests | All 5 passed, including 200 sampling seeds and official synthetic packaging |

Package SHA-256:

```text
d236c64869488827f54012a0e10de68aad05c6d97142674a561969483618371e
```

Generation and validation took approximately 11.5 minutes on the original workstation. Packaging may require roughly 40-60 GB of memory and tens of GB of temporary disk space. Large local artifacts are excluded from Git and can be rebuilt.

## Published leaderboard result

The user supplied the successful CLI output on October 3, 2026. Upload took 13 minutes 10 seconds; scoring finished with status **published**. The scores were subsequently confirmed through the CLI API; rank was 1141 at that query. The table preserves rank 1140 at publication, and rank can change as other teams submit.

| Metric | Published score |
|---|---:|
| Rank at publication | 1140 |
| Overall | **-0.276289857599391** |
| Perturbation discrimination (`pds`) | 0.019643814695154346 |
| Expression accuracy (`mse`) | 0 |
| DE log-FC accuracy (`nmae`) | 0.016005525603025445 |
| DE direction fidelity (`fid`) | **-1.6653311841128977** |
| DE direction reach (`reach`) | 0.03387419826161333 |
| DE significance overlap (`jac`) | -0.06193150004324146 |

Partition: `val`. Panel: `vcc2026-val-1`. Anchor set: `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`.

Scores are scaled against a competition baseline (0) and a split-half experimental replicate (1); the overall score is the unweighted mean of the six components. Expression accuracy is clamped to [0, 1], so its reported zero is a score floor, not zero prediction error. Direction fidelity penalizes incorrect signs and insufficient numbers of predicted responding genes. See the [official metric definitions](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics.md).

The model performed below baseline overall, principally because of direction fidelity. Poor transfer from K562 to the anonymous contexts is a plausible explanation, but the aggregate scores cannot distinguish wrong directions, insufficient DE calls or other modeling problems. The subsequent audit inspected all 900 groups and added separate experiments and a candidate zero-aware simulator. V2 was subsequently generated and published; its score and reproduction steps are in [V2_SUBMISSION.md](V2_SUBMISSION.md). The earlier audit is preserved in [MODEL_IMPROVEMENT.md](MODEL_IMPROVEMENT.md).

## Project files

| File/directory | Purpose |
|---|---|
| `build_submission.py` | Control analysis, reference transfer, integer count simulation, full validation and official packaging |
| `download_reference.py` | Resumable public reference download with pinned file ID and MD5 verification |
| `test_submission.py` | Sampling, simulation, mapping, streaming and official packaging checks |
| `requirements-competition.txt` | Versions used by the project environment |
| `manifest.json`, `gene_names.csv`, `pert_counts.csv` | Small official panel specifications needed for reproduction |
| `analysis/` | Inventory, statistics, diagnostic figure, model/validation records, tests and published score record |
| `reference/provenance.json` | Original-author reference location, checksum and CC BY 4.0 attribution |
| `reference/rejected_essential_reference.json` | Evidence that the earlier essential-gene source was not used |
| `reference/SOURCE.md` | Explanation of the retained source and excluded exploratory downloads |

Tracked analysis records include the per-gene and per-target CSVs, group statistics, source coverage and agreement diagnostics, `run_parameters.json`, official prep output, `prediction_validation.json`, `submission_ready.json` and `leaderboard_result.json`. The inventory records hashes from the original review; the README and ignore file have since changed.

The Git ignore rules exclude control/reference matrices, submission packages, local environments, extracted third-party page/paper text, downloaded website assets, and abandoned reference caches. Files remain available locally. API credentials belong in the CLI's OS keychain, not in this repository.

## Reproduce (Windows Command Prompt)

The examples below use **Command Prompt (`cmd.exe`)**, matching the successful upload session. In PowerShell, prefix executable invocations with `&` if needed and use `Set-Location` to change directory.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
```

The existing `.competition_env` is already installed on the original workstation. On another machine, use Python 3.14 to create the environment and install the recorded versions:

```cmd
python -m venv .competition_env
.\.competition_env\Scripts\python.exe -m pip install -r requirements-competition.txt
```

Obtain the official validation controls from the [competition download page](https://virtualcellchallenge.org/downloads) or the [official CLI guide](https://vcc-cli-wiki.virtualcellchallenge.org/), and place `context_A.h5ad`, `context_B.h5ad` and `context_C.h5ad` in the project root alongside the tracked panel specifications. Keep all labels and gene order intact. The saved browser pages and unrelated PDF are local review inputs; they are not needed to generate predictions from the controls.

```cmd
.\.competition_env\Scripts\python.exe test_submission.py
.\.competition_env\Scripts\python.exe download_reference.py
.\.competition_env\Scripts\python.exe build_submission.py
```

The downloader retrieves only the compact 374,587,922-byte `K562_gwps_raw_bulk_01.h5ad`, pinned to Figshare file 35774443 and MD5 `4570b53c9d62ff6df281e622f0350060`.

Additional modes:

```cmd
.\.competition_env\Scripts\python.exe build_submission.py --analysis-only
.\.competition_env\Scripts\python.exe build_submission.py --force
.\.competition_env\Scripts\python.exe build_submission.py --package-only --force
```

`--analysis-only` analyses controls without generating predictions. `--force` intentionally rebuilds/replaces existing outputs. `--package-only --force` validates and repackages the existing H5AD. The script does not upload. A rebuild rewrites generated analysis records, including the readiness record; the separate `leaderboard_result.json` preserves the first submission's historical scores.

## Submit (Windows Command Prompt)

**The first package is already submitted.** These commands document the completed workflow and future submissions. Obtain an API token from the [credentials page](https://virtualcellchallenge.org/app/credentials); paste it into the interactive login prompt.

```cmd
.\.competition_env\Scripts\python.exe -m vcc login
.\.competition_env\Scripts\python.exe -m vcc whoami
.\.competition_env\Scripts\python.exe -m vcc submit .\prediction.vcc -m "Public CRISPR transfer" --wait
```

Keep the terminal open through upload and scoring. If an upload is interrupted, resume it with:

```cmd
.\.competition_env\Scripts\python.exe -m vcc submit .\prediction.vcc --resume --wait
```

The documented allowance is two daily submissions, with one in flight. Avoid resubmitting the unchanged package merely to reproduce this result.

## Sources and limitations

- [Official competition](https://virtualcellchallenge.org/) and [CLI reference](https://vcc-cli-wiki.virtualcellchallenge.org/cli-reference).
- [Official scoring implementation and metric definitions](https://github.com/ArcInstitute/cell-eval2).
- [Original Replogle et al. processed Perturb-seq release](https://plus.figshare.com/articles/dataset/20029387), licensed [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); pinned metadata and provenance are tracked.

Coverage is partial, tissue/assay differences can invalidate transferred effects, and the original proportional simulator cannot create counts for a gene absent in a sampled input cell. The v2 zero-aware simulator addresses that limitation for genes expressed elsewhere in the context controls. Its Poisson noise approximation does not capture all biological or technical variation. File validation and passing tests establish pipeline correctness, not predictive accuracy; the published leaderboard result is the available performance evidence.

## Subsequent model audit

[MODEL_IMPROVEMENT.md](MODEL_IMPROVEMENT.md) records the step-by-step audit, measured failure mechanisms, component changes and next reference-data work. Of 900 groups, 629 call no significant downstream genes against provided controls; all contexts average below one call per target. Relaxing source suppression mostly adds changes also seen in untreated source groups. Ten focused tests now pass.

New scripts: `audit_prediction.py`, `audit_source_stages.py`, `pilot_orion.py` and `plot_model_audit.py`. The generator supports isolated `--output-dir` experiments, configurable `--transfer-shrinkage` and `--noise-z`, and optional `--upward-mode zero-aware`. Original settings remain the defaults.

A 699 MB Orion pilot was SHA-256 verified and compared on a common gene axis. Its guide library covers all targets, but the pilot has inadequate cells per target for fitting responses. A separate public CD4 T-cell summary covers 291 targets and includes independent donor/guide agreement. These are candidates for the next data step, not evidence of an improved leaderboard score.

![Model audit](analysis/model_audit/audit_overview.png)

Changed model parameters require a separate output directory, even with `--force`, to protect the first-submission files.

## Second submission candidate: donor-calibrated CD4 hybrid

The candidate is generated separately in `experiments/v2_cd4_hybrid/`. [V2_SUBMISSION.md](V2_SUBMISSION.md) records its source extraction, independent donor calibration, settings, validation and upload command. It combines eligible CD4 responses with K562 fallbacks, covering 284 targets and 11,163 measured genes across sources. Fifteen pipeline/model tests pass. It is published as `ENmII1VynYrgpB2KICPC`, with overall score **-0.223325** versus **-0.276290** for the first submission (+0.052964). Publication rank was 1105. Direction fidelity remains the largest negative component.

Reproduce source extraction and generation in Command Prompt:

```cmd
.\.competition_env\Scripts\python.exe extract_cd4_reference.py
.\.competition_env\Scripts\python.exe -m unittest test_submission test_cd4_model -v
.\.competition_env\Scripts\python.exe build_submission.py --output-dir experiments/v2_cd4_hybrid --reference-kind cd4-hybrid --transfer-shrinkage 0.7 --noise-z 1.96 --upward-mode zero-aware
.\.competition_env\Scripts\python.exe validate_v2.py
```

V2 is already published. Its completed upload command is documented below; a changed future candidate should use its own output folder:

```cmd
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v2_cd4_hybrid\prediction.vcc -m "Donor-calibrated CD4 hybrid v2" --wait
```

The generator validates and packages locally; it does not upload. Original root-level prediction files remain the first submission.

## Single CSV of supplied competition controls

`competition_controls_combined.csv` combines every supplied raw count and control-cell identifier, with genes as rows and context-prefixed cell columns. It also contains target-list flags and panel annotations. See [COMBINED_CSV.md](COMBINED_CSV.md) for the schema, training limitations and streaming loader. Regenerate with `export_competition_csv.py`; verification is recorded in `analysis/combined_csv_export.json`. The large CSV and temporary cache are excluded from Git.

## Per-cell similarity to five external references

`prepare_reference_controls.py` and `analyze_cell_reference_correlations.py` compare every A/B/C control cell with pooled non-targeting expression profiles from K562, HCT116, HEK293T, H1 and CD4. The primary Pearson comparison uses log1p(CPM) on a single common measured gene axis, with Spearman and other sensitivity analyses. All 55,200 competition cells are retained. External references are batch-limited or deterministically sampled, as detailed in [CELL_REFERENCE_CORRELATIONS.md](CELL_REFERENCE_CORRELATIONS.md).

Results, per-cell CSVs, provenance, and a comparison figure are saved under `analysis/cell_reference_correlations/`. This checks basal similarity; it does not establish anonymous cell identities or validate knockdown-response transfer.

The completed five-reference comparison uses 7,447 common genes. Median cell-to-reference correlations favor CD4 for A (0.4821) and H1 for B/C (0.4635/0.4766). All 55,200 output rows passed independent identifier, summary and checksum validation, with 45 actual-cell/reference spot checks and three focused tests. See [CELL_REFERENCE_CORRELATION_RESULTS.md](CELL_REFERENCE_CORRELATION_RESULTS.md) for the full tables and sensitivity analyses.

## Rank correlations restricted to the 300 target genes

`analyze_target_rank_correlations.py` repeats the per-cell comparison with raw values and pairwise-complete Spearman correlation, using exactly the targets in `pert_counts.csv`. No pseudocounts, transforms, normalization or missing-value imputation are applied. Measured zeros are retained as tied ranks. Each coefficient has a valid-gene count and status; K562 covers 231 targets and the other references cover 300.

See [TARGET_RANK_CORRELATIONS.md](TARGET_RANK_CORRELATIONS.md) for exact reproduction steps and [TARGET_RANK_CORRELATION_RESULTS.md](TARGET_RANK_CORRELATION_RESULTS.md) for regenerated results. Numeric outputs are in `analysis/target_rank_correlations/`.

All 55,200 rows passed validation, including 45 independent original-count comparisons and 64 pooled-profile pairs; three focused tests passed. Median target-panel correlations favor CD4 for A and H1 for B. C is close between CD4 and H1, with the median-cell and pooled summaries favoring different references. See the method document for the completed-run interpretation.

## Cell-level dispersion on the target panel

`analyze_target_cell_dispersion.py` extends the 300-target Spearman analysis with per-cell external-reference dispersion, leave-one-out correlation to the cell's own context, and exact per-cell mean/SD against every other A/B/C cell. It computes all-pair means/SDs without storing a full pair matrix; 120,000 sampled pairs provide individual examples and approximate percentiles. Five external references are additionally compared on the same 231 genes. Raw-count variance for each target, count-depth associations and guide summaries are included.

All 55,200 cells are included. Within-context cell-pair mean correlations are 0.4698/0.4949/0.4761 for A/B/C, higher than cross-context means (0.2378ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦ÃƒÂ¢Ã¢â€šÂ¬Ã…â€œ0.3168). B has the largest within-context pair SD (0.0969). Correlation varies strongly with count depth, so these descriptive results do not isolate biological heterogeneity. See [TARGET_CELL_DISPERSION.md](TARGET_CELL_DISPERSION.md) for reproduction and [TARGET_CELL_DISPERSION_RESULTS.md](TARGET_CELL_DISPERSION_RESULTS.md) for tables and the figure. Outputs are in `analysis/target_cell_dispersion/`.

Validation passed for all 55,200 cell rows, 120,000 sampled pair rows and 900 gene-variance rows. Independent checks included 27 explicitly enumerated all-peer mean/SD comparisons, 45 shared-axis reference correlations, nine leave-one-out correlations and 36 individual pairs against SciPy. Five focused tests passed; evidence is in `analysis/target_cell_dispersion/validation.json` and `test_results.txt`.

## V3 context-specific response ensemble

[V3_SUBMISSION.md](V3_SUBMISSION.md) records the new strategy and exact reproduction/submission commands. It combines uncertainty-calibrated CD4 condition responses, direct public 2025 H1 responses for 13 targets, and conservative K562 fallbacks. Basal expression-fit weights differ between A/B/C; they are not inferred cell-type fractions. Weights are renormalized per target and measured gene, and uncertain CD4 condition weights are pulled halfway toward equal consensus.

Five nested target-grouped donor folds selected the same CD4 uncertainty estimator: relative source error 0.984254 versus V2 0.986857. H1 held-target/disjoint-batch selection error is 0.802987 relative to zero response. These are source-domain results, not hidden VCC accuracy estimates. All contexts have downstream estimates for 284 targets. Each context uses a common set of 400 distinct, guide/depth-balanced controls across its targets to reduce baseline sampling differences. Generation uses the existing zero-aware integer simulator and is isolated in `experiments/v3_context_ensemble/`.

Root code: `extract_h1_responses.py`, `fit_context_ensemble.py`, `test_context_ensemble.py`, `validate_v3.py`, `audit_v3_effects.py` and `plot_v3_strategy.py`, plus the ensemble integration in `build_submission.py`. Source/model axes, fold separation, checksums and control templates passed validation; all 19 generation/model regression tests passed. The 3,305,861,120-byte package passed official preparation, complete count/label validation and 333,594 independent gene-count expectation checks across 18 selected groups. It is ready at `experiments/v3_context_ensemble/prediction.vcc`; package checksum and validation evidence are recorded in the candidate's `analysis/` directory. No new leaderboard score is available until upload and scoring.

## Complete-matrix rank clustering of A/B/C controls

`cluster_full_rank_matrix.py` clusters all 55,200 original cells together on all 18,533 measured genes using raw average-rank Spearman geometry. Zeros remain valid tied measurements; there are no pseudocounts, log transforms, gene filters or count normalization. Context labels are withheld from fitting. Full-gene spherical clustering is evaluated at K=2..6 and K=3 is repeated with a second seed. PCA is used only to display the fitted cells.

Both independent three-cluster fits recover A, B and C perfectly: each cluster contains all 18,400 cells of one context, with no mixing (ARI/NMI/purity 1.0). All-peer mean correlations are approximately 0.63 within contexts and 0.44ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Å“0.51 across contexts. Higher K mostly isolates single cells; K=6 additionally separates 1,013 B cells. The 1,500-cell silhouette sample omits singleton clusters, so its tied K=3/4/5 scores do not prove that the full partitions are equally good. Count-depth gradients remain within the clusters, and context separation may include batch effects.

[Methods and reproduction](FULL_RANK_CLUSTERING.md) and [results](FULL_RANK_CLUSTERING_RESULTS.md) document every step. `analysis/full_rank_clustering/per_cell_full_rank_clusters.csv` contains the original cell IDs, every fitted cluster assignment, depth/detection, centroid correlations/margins, exact all-peer context means and eight display PCs. Other outputs include complete composition tables, convergence histories, depth/guide diagnostics, source/output hashes, saved centroids and a PNG/PDF figure.

```cmd
.\.competition_env\Scripts\python.exe -m unittest test_full_rank_clustering -v
.\.competition_env\Scripts\python.exe cluster_full_rank_matrix.py --seed 20261004 --threads 12
.\.competition_env\Scripts\python.exe validate_full_rank_clustering.py
```

Three focused tests and an independent validator cover original identifiers, gene order, SHA-256 provenance, metric/composition consistency, raw-count reranking, saved-model assignments, PCA coordinates and full-gene sample silhouettes. These results characterize the supplied controls; they do not establish anonymous cell identities or perturbation-response accuracy.

## External controls mapped into the A/B/C rank space

`extract_external_mapping_controls.py` recovers the previously selected 4,393 external observations from pinned local files and checksum-verified cached CD4 ranges, without further downloads. These comprise 503 K562 pseudobulk groups, 891 HCT116 cells, 1,079 HEK293T cells, 1,536 H1 cells and 384 CD4 pseudobulk groups. Source-unit sums recover the previously verified pooled control profiles.

A valid shared space uses **7,447 jointly measured genes**, because the external releases do not all measure the original 18,533-gene axis. `map_external_rank_controls.py` reranks all 55,200 competition cells on that common axis, fits label-blind spherical clusters and eight PCs using competition controls only, freezes the model, and projects the external observations. Zeros are retained; there is no pseudocount, log, CPM or missing-symbol imputation. The shared-axis clusters still recover A/B/C perfectly (context ARI and repeat ARI 1.0).

Closest-context fractions are CD4ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢A 100%, H1ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢C 98.96%, K562ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢A 100%, HCT116ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢C 97.64% and HEK293TÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢A 97.03%. K562/CD4 fractions count pseudobulk observations, not individual represented cells. All external observations fall outside every context's empirical 95th-percentile region in eight PCs, so these are relative basal matches rather than literal membership or established cell identities. The stronger pooled-profile matches are CD4ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Å“A (Spearman 0.717) and H1ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Å“C (0.702); H1ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Å“B remains 0.656. Other source matches are weaker.

[EXTERNAL_RANK_MAPPING.md](EXTERNAL_RANK_MAPPING.md) documents sources, coverage, unit definitions, frozen projection, distribution-shift checks, output schemas and reproduction. [EXTERNAL_RANK_MAPPING_RESULTS.md](EXTERNAL_RANK_MAPPING_RESULTS.md) records all fractions, pooled comparisons, CD4 condition sensitivity and figures. Complete observation-level coordinates, metadata and similarities are in `analysis/external_rank_mapping/external_control_units_mapped.csv`; all 55,200 competition coordinates are also saved. The original full-gene clustering and prediction files are retained.

```cmd
.\.competition_env\Scripts\python.exe -m unittest test_external_rank_mapping -v
.\.competition_env\Scripts\python.exe extract_external_mapping_controls.py
.\.competition_env\Scripts\python.exe map_external_rank_controls.py --seed 20261004 --threads 12
.\.competition_env\Scripts\python.exe validate_external_rank_mapping.py
```

The independent validator checks all identifiers, source rows, hashes and assignment summaries, reranks original competition and external observations, verifies frozen model projections and shift distances, and recalculates all pooled-profile Spearman correlations. Focused tests check measured-zero/duplicate-symbol handling, fixed projection and corrupt/missing offline ranges.


## Verified 300-target ZIP, target-only clustering and lineage markers

The supplied `vcc_2026_300_target_matrix.zip` contains all 300 target genes and 18,400 controls per A/B/C context. `analyze_target_archive.py` verifies all **16,560,000 counts exactly against the original H5AD matrices**, with zero mismatches; ZIP CRC, gene coverage/order and integer counts also pass. The CSV has no cell IDs; IDs are attached from the completely verified within-context source order. The original ZIP is kept locally as an input with its checksum recorded in the audit.

Blind clustering of all 55,200 controls uses average tied ranks on only these 300 raw-count genes, keeping zeros and adding no pseudocount or normalization. Both three-cluster seeds give the same partition: **55,190 cells (99.9819%) agree with A/B/C**, with ten low-depth B cells joining A (context ARI 0.999457, repeat ARI 1.0). Those cells have 737ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Å“1,298 full-library counts and 8ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Å“27 detected target genes; the earlier full-gene clustering correctly groups them with B.

Full-library CPM and pairwise rank AUC identify 37 strong dataset-relative enriched genes under explicit thresholds: 20 in A, 5 in B and 12 in C. Leading markers are A: SMARCA5/PIN1/HDAC1/MARCKSL1/EZH2; B: VIM/PTMS/ACLY/ABCD1/SND1; C: S100A11/PRNP/FXYD5/RAB5A/PATL1. KIF21B in A is the only target meeting the stricter detection-restriction rule. Complete 900-row tables are saved for both fitted clusters and original contexts. Enriched expression does not imply universal tissue specificity.

An additional 50-gene original-matrix check supports a T-cell-like program in A (CD3D/CD3E/TRAC), mixed epithelial/mesenchymal-associated expression in B (VIM plus KRT8/KRT18), and an epithelial program in C (EPCAM/CDH1/keratins with very low VIM). S100A11 and FXYD5 alone are not diagnostic epithelial markers, and chromatin/mitosis expression alone cannot establish tissue of origin. Lineage hypotheses are documented with primary-source links in [TARGET_LINEAGE_INTERPRETATION.md](TARGET_LINEAGE_INTERPRETATION.md).

Root scripts: `analyze_target_archive.py`, `test_target_archive.py`, `validate_target_archive.py`. [TARGET_ARCHIVE_AUDIT.md](TARGET_ARCHIVE_AUDIT.md) describes every step, inputs, marker rules, limitations and output schemas; [TARGET_ARCHIVE_RESULTS.md](TARGET_ARCHIVE_RESULTS.md) contains numeric results and figures. Per-cell IDs, all marker metrics, count-depth sensitivity, models, audit hashes and PNG/PDF figures are under `analysis/target_archive_audit/`.

```cmd
.\.competition_env\Scripts\python.exe -m unittest test_target_archive -v
.\.competition_env\Scripts\python.exe analyze_target_archive.py --seed 20261004 --threads 12
.\.competition_env\Scripts\python.exe validate_target_archive.py
```

Three focused tests passed. Independent validation repeats every one of the 16,560,000 original count comparisons and verifies all 55,200 cell identities/depths, all 1,800 cluster/context marker rows, all 150 broader marker summaries, 70 independently reranked/model-checked cells and all saved output/script checksums. Evidence is in `analysis/target_archive_audit/validation.json`. This analysis characterizes the provided controls; it does not generate or upload a new submission.

## Neural networks for the 300-gene panel

Compared shallow, deep, bottleneck and residual neural configurations separately in A/B/C with normalized-expression and direct raw-count objectives (24 candidates). Validation guide groups selected architectures and stopping epochs; frozen models were compared on the existing ridge test set.

Selected neural raw-count MSE reductions against the depth-only baseline: **A: 31.94%**, **B: 25.97%**, **C: 28.79%**. Selected normalized-expression reductions against the mean profile: **A: 7.99%**, **B: 10.87%**, **C: 8.38%**.

Saved separate normalized/raw models, both held-out prediction types and full raw reconstructions for all 55,200 controls. Five focused tests, six standalone inference checks and independent artifact/matrix checks passed. Methods and commands: [GENE300_NEURAL_MODELS.md](GENE300_NEURAL_MODELS.md). Architecture comparison and figures: [GENE300_NEURAL_RESULTS.md](GENE300_NEURAL_RESULTS.md). Root scripts include `train_gene300_neural.py`, `predict_gene300_neural.py`, `test_gene300_neural.py`, `validate_gene300_neural.py` and `document_gene300_neural.py`.

## Separate models from the 300-gene panel

Trained separate A/B/C models using 300 measured genes to predict the other 18,233 genes. Held out complete control-guide groups for validation and test; saved trained weights, predictions for all 55,200 controls, 8,400 held-out predictions and detailed gene/cell metrics.

Held-out normalized-expression MSE reductions against each context's training mean: **A: 7.20%**, **B: 9.70%**, **C: 7.40%**. Raw-count conversion is compared separately with a depth-only baseline; see the results report for its limitations. Full-data reconstructions are in-sample. These observational control models are separate from competition submission artifacts.

Methods and commands: [GENE300_MODELS.md](GENE300_MODELS.md). Measured results and figures: [GENE300_MODEL_RESULTS.md](GENE300_MODEL_RESULTS.md). Root scripts: `train_gene300_models.py`, `predict_gene300_models.py`, `validate_gene300_models.py`, `test_gene300_models.py`, `smoke_gene300_inference.py` and `document_gene300_models.py`.

## V6 fidelity-guarded submission

Created and independently verified `experiments/v6_fidelity_calibrated/prediction.vcc` (3,309,619,200 bytes), with all 900 groups, 360,000 cells and 18,533 genes. Frozen candidate: `blend25_v4`, settings `{"blend": 0.25, "calibration": 0.0, "dual": false}`. A retains V4 effects; B/C updates are bounded and cannot reverse V4 signs or create previously zero effects. Neural templates were evaluated and withheld after distribution/perturbation-evidence checks.

Seven tests, disjoint source/control identity checks, all-coordinate model replay, the full count/library scan, every packaged count/index/pointer, all source-row provenance checks and 12 exact full-group replays passed. Official VCC packaging passed. [V6_SUBMISSION.md](V6_SUBMISSION.md) records actual source benchmarks, limitations, settings, hashes and Command Prompt reproduction/upload commands. Published as `6MCjBWpPoCt5Sbw5CMu3`, overall **-0.056790** versus V4 **-0.050674**. Fidelity and reach decline despite gains in PDS, LFC accuracy and overlap. [V6_LEADERBOARD_RESULTS.md](V6_LEADERBOARD_RESULTS.md) records the API-confirmed raw/scaled metrics. V4 was best at that publication.

## Next-submission strategy after V5 and neural modeling

That earlier V6 strategy used V4 as its starting point. Compare response updates and count generation separately before adding a neural cell-state template. Across all 900 frozen effect groups, V5 changes only 22 A coordinates; fewer than 1% of sizeable effects in both versions reverse sign in B/C. This is descriptive model agreement, not correctness against hidden responses. [V6_MODEL_STRATEGY.md](V6_MODEL_STRATEGY.md) specifies the ablations, source/guide isolation, distribution and direction-yield checks, neural limitations and adoption rules. Root diagnostic: `analyze_v6_strategy.py`; records: `analysis/v6_strategy/`. That proposal and diagnostic preceded the completed V6 experiment and publication documented above.

## V5 atlas-calibrated submission

A new candidate is generated and verified at `experiments/v5_atlas_calibrated/prediction.vcc` (3,275,069,440 bytes; 3.05 GiB). It covers all 900 groups, 360,000 cells and the exact 18,533 genes. Official VCC packaging and independent source/model/count verification passed. Published as `BPI0n52HTsG0EYTFOf9H` with overall **-0.113637**; V4 was best at that publication, at **-0.050674**.

V5 expands public H1 coverage to 25 source targets (25 QC-accepted direct responses), estimates mean/pooled effects separately, generates from individual control-cell templates while preserving all original cell depths, and applies 22 source-supported promoter-neighbor corrections only in A. CD4 across-target subtraction remains zero. Five tests and 24 numerical preflight groups passed. [V5_SUBMISSION.md](V5_SUBMISSION.md) gives the upload command and hash; [V5_MODEL_STRATEGY.md](V5_MODEL_STRATEGY.md) gives every step, selected settings, verification and limitations.

Root scripts are `extract_v5_h1.py`, `calibrate_v5_promoters.py`, `fit_v5_submission.py`, `build_v5_submission.py`, `test_v5_submission.py`, `validate_v5.py` and `document_v5_submission.py`. Small evidence is tracked in the candidate analysis folder; source counts/caches and large predictions are regenerated locally. At that stage, V4 was the best of six published project submissions. V5 improves five scaled components, but direction fidelity falls sharply. [V5_LEADERBOARD_RESULTS.md](V5_LEADERBOARD_RESULTS.md) records all scores, raw metrics and the arithmetic decomposition.

## Upstream repository review and validated improvement experiments

Both cloned repositories were reviewed through their reachable commit history, historical methods and current code/tutorials. [UPSTREAM_REPOSITORY_REVIEW.md](UPSTREAM_REPOSITORY_REVIEW.md) records pinned commits, source differences, tests, scientific figures and exact reproduction commands. AtlasShift has 14 commits but only two executable-Python-changing commits; Arc Atlas has 59 main commits (62 across refs) and supplies data/access tutorials rather than a trained predictor.

Public 2025 H1 metadata extends potential direct-response coverage from 13 to **25** of our 300 targets, with 12 additional targets. All 38,176 control identifiers recur across the released splits, so validation must partition controls by identifier. Response counts/knockdown still require extraction and QC. Paired neutral-control tests show that separate mean-CPM/pooled-profile count matching reduces profile drift while one-cell templates retain variance near the controls; four-cell pooling reduces variance to roughly 30-38%. Five held-target CD4 folds select zero across-target subtraction. Independently reconstructed promoter pairs cover 74 targets, with descriptive neighboring-gene suppression evidence that requires held-out calibration.

Root scripts: `review_upstream_repositories.py`, `probe_atlas_sources.py`, `evaluate_atlas_transfer.py`, `test_upstream_review.py`, `validate_upstream_review.py`, `document_upstream_review.py`. Four focused tests and the independent validator passed; evidence is under `analysis/upstream_review/`. The next candidate should first evaluate additional H1 responses and one-cell dual-moment calibration. This review preserves V4 and does not generate or upload a new prediction package.

## V4 lineage-guided response coverage candidate

The completed candidate is ready at `experiments/v4_lineage_response/prediction.vcc` (3,309,240,320 bytes; 3.08 GiB). Official preparation and independent package verification passed for all 900 groups, 360,000 cells and the exact 18,533-gene axis. Published as `w5criCkgkCWj1QOgP08E` with overall score **-0.050674**. Rank was 977 at publication and 975 at the recorded API check; ranks change over time. [V4_SUBMISSION.md](V4_SUBMISSION.md) records the package hash, completed checks, actual generated-response audit and Command Prompt upload command.

The model uses the marker-supported T-cell-like program to prefer CD4 responses for A, strengthens donor-supported target effects to address response undercoverage, and uses public H1 direct responses plus held-target-calibrated transfer/shared responses for B/C. It applies explicit basal gene/target expression gates and independently draws 400 real controls per target with guide/depth balance. Source selection improves direction-coverage proxies with a documented squared-error tradeoff; it does not establish a better hidden leaderboard score. [V4_MODEL_STRATEGY.md](V4_MODEL_STRATEGY.md) documents every step and exact reproduction.

New root code is `fit_lineage_submission.py`, `test_lineage_submission.py`, `validate_v4.py` and `audit_v4_strategy.py`; `build_submission.py` adds the isolated `lineage-response` mode and complete original-row provenance. All 23 tests passed. Independent validation recalculates all public-source response components, reassembles all 16,679,700 model effects, verifies all 360,000 original control-row/guide/depth records and checks 389,193 gene-count expectations with exact selected-cell generation replay. Model/package provenance and evidence are retained in the candidate analysis folder.

## Published validation result

The first, V2 and V4 scores use the same validation panel and anchor set. No V3 submission score is recorded, so there is no V3 leaderboard comparison.

| Scaled component | First | V2 | V4 | Change from V2 |
|---|---:|---:|---:|---:|
| overall | -0.276290 | -0.223325 | -0.050674 | +0.172652 |
| pds | 0.019644 | 0.103042 | 0.255302 | +0.152260 |
| mse | 0.000000 | 0.000000 | 0.000000 | +0.000000 |
| nmae | 0.016006 | 0.032352 | 0.067539 | +0.035186 |
| fid | -1.665331 | -1.470413 | -0.683855 | +0.786558 |
| reach | 0.033874 | 0.041380 | 0.086849 | +0.045469 |
| jac | -0.061932 | -0.046314 | -0.029878 | +0.016436 |

Overall improved by **+0.172652** versus V2 and **+0.225616** versus the first entry. All five non-expression scaled components improve. Direction fidelity accounts arithmetically for 75.9% of the overall gain; this is a score contribution, not causal attribution to a model change.

Raw expression error worsens from 1.014294 to 1.247733 (+23.0%), although the scaled expression score stays at zero. Zero is a clamped leaderboard score, not zero prediction error. Raw direction-fidelity yield rises from 0.074675 to 0.309235; it is response-budget-normalized yield rather than sign precision. Raw significant-set Jaccard rises from 0.013545 to 0.019731, leaving substantial response-set mismatch.

At that stage, V4 was the best of the six published project submissions on this validation panel; the later V5 result is documented above. At V4 publication, overall was below the zero-score anchor, with direction fidelity the largest negative component. The package changes lineage policy, effect strength, transfer priors and control sampling together; aggregate scores cannot isolate which change helped, which context remains weak, or whether residual fidelity loss comes from missed responses or wrong directions. The next experiments should balance response coverage with expression error using source validation and controlled ablations.

Records: `leaderboard_result.json`, `published_api_status.json`, `leaderboard_comparison.csv` and `raw_metric_comparison.csv` in the candidate analysis folder. `submission_ready.json` preserves the generation-time snapshot; the separate publication record records completed upload and scoring. Upload took 13 minutes 20 seconds. [Official scoring definitions](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics.md).

Record or refresh the read-only API confirmation and regenerate these tables with the root script `record_v4_result.py`:

```cmd
.\.competition_env\Scripts\python.exe record_v4_result.py
```

The supplied publication rank is preserved separately from the current API rank. This command does not upload, regenerate predictions, or alter the package.


## Next performance experiments (2026-10-07)

[PERFORMANCE_IMPROVEMENT_PLAN.md](PERFORMANCE_IMPROVEMENT_PLAN.md) compares our implemented response models with current public AtlasShift, Sidechain and STATE methods. The proposed sequence prioritizes full measured CRISPRi source coverage, adaptive response calibration, nested response-head tuning and isolated count-emission experiments. This section records the October 7 research plan. The subsequent V10 implementation and upload-ready artifact are documented below. V4 was the best published candidate when this October 7 plan was written.

<!-- V10_VERIFIED_SUBMISSION -->

## V10 published result

V10 is published as `CYpAqJP7bUYqjA2ouD7U`, overall **0.116754121**, an improvement of **+0.167428032** over V4. Official prep and independent full-payload validation passed. The pooled dose-corrected response and compact emitter improve discrimination and fidelity, while reach and fold-change accuracy decline. Compact neutral-control bias remains. See [V10_LEADERBOARD_RESULTS.md](V10_LEADERBOARD_RESULTS.md), [V10_STRATEGY.md](V10_STRATEGY.md) and [V10_SUBMISSION.md](V10_SUBMISSION.md).

<!-- V11_VERIFIED_SUBMISSION -->

## V11 published result

V11 is published as `bmIAGveuEvFiYQTvz0Ma`, overall **0.137003974**, rank 529 at publication, improving on V10 by **+0.020249853**. Five of six scaled components improve; expression accuracy declines. Official preparation and independent full-payload/replay/provenance audits passed. Source consensus rescue and higher dispersion improve coverage together; neutral-control bias remains. See [V11_LEADERBOARD_RESULTS.md](V11_LEADERBOARD_RESULTS.md), [V11_STRATEGY.md](V11_STRATEGY.md) and [V11_SUBMISSION.md](V11_SUBMISSION.md).
