<!-- project-progress-20261023:begin -->
## Current status: final V91 uploading (23 October 2026)

The user reports that the validated D/E/F submission is uploading. Completion and entry ID have not yet been received; final scores remain hidden. **V90 remains our best published validation score: 0.15214377976638518.**

[Complete project progress and all 22 validation results](PROJECT_PROGRESS.md) covers the data analyses, models, rejected experiments, final source reconstruction, validation and reproducibility. [Upload status snapshot](analysis/project_progress_20261023/upload_status.json) records the user report separately from frozen pre-upload receipts. Older readiness/status blocks below retain their historical meaning. Do not regenerate the active uploading file.
<!-- project-progress-20261023:end -->

# Final-round preparation

The final bundle was inspected on October 22, 2026. It has 18,400 non-targeting controls in each of **D, E and F**, 300 perturbation targets with 400 predicted cells each, and the same 18,533-gene axis as validation. Only **129 targets overlap** the old list; **171 are new**. Keep labels exactly as provided. Validation predictions must not be relabelled D/E/F.

The project-local VCC CLI was upgraded from 0.2.2 to **0.3.0**. The final `manifest.json` must remain beside `final/pert_counts.csv`. Packaging must use both final CSVs, with D/E/F explicitly checked. The official guide confirms that the phase is blind and only the latest final submission counts: https://vcc-cli-wiki.virtualcellchallenge.org/ . The deadline supplied by the organizer is **November 5, 2026 at 23:59 UTC** (November 6 at 00:59 Europe/Oslo).

## Reproduction and current stages

Run from the project root with `.competition_env\Scripts\python.exe -X utf8`:

1. `inspect_final_inputs.py`: validates every final control count, original gene/context labels and library/guide structure; hashes all input files and writes `experiments/final_v91/analysis/input_audit.json` plus final-control moments.
2. `record_final_transition.py`: stores the user-supplied published V90 result separately from frozen pre-upload evidence.
3. `collect_final_cd4.py`: extracts only the 171 new targets from the same ETag-pinned merged and disjoint-donor DE summary files; retains original validation-target rows exactly. Source identity, byte bounds and checksums remain enforced.
4. `collect_final_orion.py HCT116` and `collect_final_orion.py HEK293T`: extend the same pinned full-release raw-cell source with only new targets. The original reviewed sampler, missing-fragment offset repair and complete-source library checks are reused in isolation; old source caches are not overwritten.

5. `collect_final_h1.py`: extracts newly required targets from the three public 2025 H1 splits using pinned byte-range downloads. The same verified H1 NTC rows are retained once when summarizing; they are not duplicated across public splits.
6. `audit_final_replay.py A`, `B`, `C` and `test_final_count_replay.py`: replay archived response/count arrays under the working runtime. **All nine V24/V66/V90 response arrays are bitwise identical; three full 400-cell archived populations also regenerate exactly.** These technical checks do not establish final biological accuracy.
7. `test_final_cli.py`: offline round-trip of D/E/F manifest inference, panel stamp and perturbation-list fingerprint under CLI 0.3.0. The toy file is temporary, deleted and never uploaded.
8. `run_final_pipeline.py`: waits for completed collection/replay receipts; rebuilds and audits the final source cache; freezes the exact method and code/input hashes; fits real and shuffled source priors; fits D/E/F responses; audits numerical/identity guards; runs 12-group count preflight; streams all 360,000 predictions; packages with official CLI 0.3.0 and independently checks the complete payload and twelve regenerated groups. It never uploads or pushes.

The new public response sources must be audited before model fitting. Unsupported target/output coordinates must abstain rather than become fabricated source measurements. New final-control cells supply basal state and templates; no final treatment labels are available. No claim of final score improvement is possible from control checks.

The old temporary Python3.12 scientific environment is currently incomplete. This preparation uses the working project-local Python3.14 runtime, with versions recorded in the input audit. Validation replay is required before releasing adapted response/count predictions; exact regeneration must not be assumed across runtimes.

**The final upload file is ready at `experiments/final_v91/prediction.vcc`. Official prep and independent full-payload validation passed. No final submission has been uploaded.** See [FINAL_V91_RESULTS.md](FINAL_V91_RESULTS.md) for the checksum and CMD submission command. Large controls, extracted sources and generated arrays remain local and ignored by Git. Root scripts and compact audit evidence are versioned.

## Completed source work and basal diagnostics

The full-release Orion extension contains **30,795 HCT116 cells across 109 batches** and **42,155 HEK293T cells across 223 batches** (including controls), each covering 169 of the 171 new target names. All source-library totals and repaired offsets were checked. Coverage after the 30-cell response-summary requirement is reported separately in `target_coverage.csv` after assembly; merely having a target name in a source is not sufficient evidence.

Pairwise finite-value Spearman correlations of final control mean CPM, without pseudocounts, give nearest available external controls: D–CD4 Stim48hr **0.850588**, E–H1 **0.776744**, F–H1 **0.802042**. These are basal similarities, **not cell-line identities or validation of intervention-response transfer**. See `basal_rank_correlations.csv` and `context_diagnostic.json`.

The source-only final model retains the exact published method choices (V24 `cdf50_de`, V66 `assay75`, V90 `interval25`, residual10 count emitter). It refits the CD4 teacher, source centering, rank axes and finite-update covariance prior for the actual final panel. Missing source observations stay masked. Final treatment outcomes remain unavailable; public validation rank cannot forecast the final rank.

## Frozen final response source audit

The final cache contains **1,487** public target/context response rows: **638 reused bitwise** and **849 newly extracted**. It covers **298/300 final targets**. **ANG and DPY30 have no eligible source response** and keep a control-based downstream baseline with the same preset 90% on-target suppression assumption; all 400 predicted cells per context are still required for these targets. No source response is invented for missing evidence.

Source-family row counts are K562 260, CD4 625 (three correlated conditions), H1 28, HCT116 282 and HEK293T 292. Family coverage per final target is: 19 targets with five families, 180 with four, 82 with three, 14 with two, 3 with one and 2 with none. These are measurement-coverage counts, not performance estimates. Source masks and own-coordinate exclusions passed independent checks.

The added H1 evidence comprises **8,463 treated cells** for 12 newly required targets across the public 2025 train/validation/test splits. Byte-range identities and checksum checks were retained. Parallel verified-page loading accelerated shared-drive reads without changing sample selection or counts.

The scientific method was frozen before final response inference in `experiments/final_v91/analysis/plan.json`. Both original finite-update source priors were fitted and saved, including the shuffled diagnostic. Full diagnostic response outputs are retained for inspection; only the preselected V90 interval25 output is used for cells.

Completed: `summarize_final_generation.py` reproduces both preflight and all-900-group moment diagnostics. `report_final_results.py` writes the readiness report only after both official and independent receipts pass. `audit_final_repository.py` checks frozen code/index hashes, final input integrity and preservation of earlier protocols.
