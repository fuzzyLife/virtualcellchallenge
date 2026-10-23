# Project progress through 23 October 2026

This is the current navigation and evidence summary. **Final V91 is uploading according to the user; completion and entry ID have not been received.** Local preparation and validation are complete. Final scores are hidden. The best published validation result remains **V90: 0.15214377976638518**.

Earlier README sections and experiment reports are dated snapshots. Their historical statements such as "best published" or "not uploaded" describe what was known then. This document consolidates the later receipts without rewriting frozen scientific evidence.

## 1. Competition task and data understanding

The challenge provides non-targeting control cells and an intervention list, not matched treatment labels for the anonymous challenge contexts. We must predict raw integer expression counts after CRISPRi, on the supplied gene axis, with the required number of cells for every target/context. Control-expression reconstruction and causal intervention prediction are different tasks.

The validation panel uses A/B/C. The final panel uses **D/E/F**, each with **18,400 non-targeting control cells**, 46 guides with 400 cells each. The gene axis has **18,533 genes**. Each final context requires **300 targets x 400 predicted cells**, so the output has **360,000 rows, 18,533 columns and 900 target/context groups**. Cells are the rows; genes are the columns. `pert_counts.csv` specifies targets and required prediction-cell counts; it is not an observed perturbed expression matrix.

The six metrics are reference-scaled within context. Zero corresponds to the context-mean baseline; one corresponds to the experiment replicate reference. The overall is their mean across contexts and metrics. CLI MSE/NMAE score columns are scaled accuracies, with higher favorable; raw error columns have different units. The perturbed gene itself is excluded from downstream scoring, so suppressing only the target cannot solve the task.

We inspected the official evaluation and FAQ, checked format/size/library requirements and repaired local scorer issues, including discrimination origins, count-depth adaptation, eligibility, tied ranks and significance arithmetic. Public tests do not reproduce private challenge anchors. See [task analysis](COMPETITION_ANALYSIS.md), [evaluation audit](EVALUATION_CRITERIA_AUDIT.md), [FAQ audit](FAQ_AUDIT.md) and [transfer-validation audit](TRANSFER_VALIDATION_AUDIT.md).

## 2. Data exploration, gene panels and observational models

- Built combined control-data exports with explicit context/gene metadata: [COMBINED_CSV.md](COMBINED_CSV.md).
- Compared control profiles with K562, HCT116, HEK293T, H1 and CD4 references; then restricted to the 300 targets and used pairwise valid-value Spearman ranks without pseudocounts: [reference correlations](CELL_REFERENCE_CORRELATION_RESULTS.md), [target-rank analysis](TARGET_RANK_CORRELATION_RESULTS.md).
- Measured per-cell correlations and within-context dispersion, clustered full expression matrices and mapped external controls into the same rank space: [cell dispersion](TARGET_CELL_DISPERSION_RESULTS.md), [full clustering](FULL_RANK_CLUSTERING_RESULTS.md), [external mapping](EXTERNAL_RANK_MAPPING_RESULTS.md).
- Audited the supplied 300-gene ZIP against original values and gene coverage, tested its context separation and examined cluster markers: [archive audit](TARGET_ARCHIVE_RESULTS.md), [lineage interpretation](TARGET_LINEAGE_INTERPRETATION.md). VIM or epithelial-associated markers alone do not uniquely identify the anonymous cell lines.
- Trained separate per-context 300-input-gene reconstruction models, including ridge and several neural structures: [linear results](GENE300_MODEL_RESULTS.md), [neural results](GENE300_NEURAL_RESULTS.md). These measured observational coexpression; modifying an input gene is not validated CRISPRi response inference.
- Produced control/prediction cell projections and examined population collapse/dispersion: [V16 projection](V16_CELL_PROJECTION.md). A broader scatter cloud alone is not a better intervention prediction.

## 3. All recorded validation submissions

There are **22 recorded published validation entries**, including the V82 recovery of identical V66 predictions. Version numbers also cover many offline experiments; they do not count submissions or independent model advances. All listed results use `val / vcc2026-val-1` and the same saved validation anchors. Compare scores rather than ranks from different dates. Tiny gains have no established statistical significance.

| Entry | Main released change | Overall | Change from previous project best | New best? |
|---|---|---:|---:|:---:|
| V1 | Public CRISPR transfer | -0.276289858 | initial | yes |
| V2 | Donor-calibrated CD4 hybrid | -0.223325462 | +0.052964396 | yes |
| V4 | Lineage-guided response coverage | -0.050673910 | +0.172651551 | yes |
| V5 | Atlas-calibrated dual moments | -0.113636731 | -0.062962821 |  |
| V6 | Fidelity-guarded lineage response | -0.056790139 | -0.006116229 |  |
| V7 | V4 fallback resampling | -0.051386751 | -0.000712840 |  |
| V9 | Activation-preserving profile calibration | -0.070447281 | -0.019773371 |  |
| V10 | Pooled intervention population transport | 0.116754121 | +0.167428032 | yes |
| V11 | Consensus coverage and dispersion | 0.137003974 | +0.020249853 | yes |
| V13 | Assay-aware knockdown transport | 0.137883134 | +0.000879160 | yes |
| V16 | Response residual dispersion | 0.141613857 | +0.003730723 | yes |
| V20 | Stack-calibrated response ensemble | 0.139011567 | -0.002602290 |  |
| V24 | Context-matched donor response | 0.142997855 | +0.001383998 | yes |
| V36 | Retrieval-contrast correction | 0.141845482 | -0.001152374 |  |
| V50 | Empirical response prior | 0.141522622 | -0.001475234 |  |
| V55 | Fixed V24/V36 half-blend | 0.138767074 | -0.004230781 |  |
| V61 | Pooled state-response transfer | 0.116763474 | -0.026234381 |  |
| V63 | Identity-guarded assay response | 0.150422795 | +0.007424940 | yes |
| V66 | Assay-strength response | 0.150882188 | +0.000459393 | yes |
| V69 | Probability-coupled count mixture | 0.095262845 | -0.055619343 |  |
| V82 | V66 recovery; identical predictions | 0.150882188 | same predictions as V66 |  |
| V90 | Conservative interval projection | 0.152143780 | +0.001261592 | yes |

Exact entry IDs, component scores, source-receipt paths and full-precision values are in [the consolidated JSON](analysis/project_progress_20261023/published_history.json) and [CSV](analysis/project_progress_20261023/published_history.csv). V69 FID/REACH were not reported and stay missing. V55 uses its subsequently saved authenticated receipt, not just the initially rounded webpage paste.

The main observed advances were stronger source coverage/lineage policy (V4), pooled intervention responses plus population transport (V10), consensus coverage/dispersion (V11), empirical control residuals (V16), context-matched CD4 responses (V24), and assay-response identity guards (V63). V66 and V90 added much smaller gains. These are bundled changes and selected comparisons, not randomized causal component estimates.

V10 improved over V4 by 0.167428032. Its fidelity and discrimination score changes account arithmetically for about 97.4% of the net gain. V20 stacking, V36 retrieval correction, V50 priors and V55 blending did not beat V24. V61 regressed. V69 changed count populations and lost 0.055619343 overall relative to V66 despite valid formatting; V82 reproduced the original V66 score with identical predictions. V90 beat V66 by only 0.001261592, improving five components while expression accuracy slipped.

Read the [early retrospective](PROJECT_RETROSPECTIVE.md), [method/performance association analysis](V54_METHOD_PERFORMANCE.md), [V65 history](V65_SUBMISSION_HISTORY.md), [V66 receipt](V66_LEADERBOARD_RESULTS.md), [V69 regression](V69_LEADERBOARD_RESULTS.md), [V82 recovery](V82_LEADERBOARD_RESULTS.md) and [V90 receipt](V90_LEADERBOARD_RESULTS.md) for the original interpretations and limitations.

## 4. Research and offline ideas, including negative results

| Work | What was tested or reviewed | Outcome and evidence |
|---|---|---|
| V3-V9 | Basal deconvolution, lineage transfer, atlas/moment calibration, neural transfer heads and resampling | Expression similarity did not establish causal weights. V7 released V4 effects with new sampling; several changes lost fidelity. [Upstream review](UPSTREAM_REPOSITORY_REVIEW.md), [V7 strategy](V7_MODEL_STRATEGY.md). |
| V10-V13 | Pooling, compact population transport, distribution variants, assay harmonization, low-rank response programs and knockdown calibration | Pooling/coverage improved hosted scores. Many distribution/program variants failed; V13 selected knockdown handling rather than a learned program network. [V12 research](V12_RESEARCH.md), [V13 strategy](V13_STRATEGY.md). |
| V14-V24 | Specialists, bounded ensembles, residual populations, direction/margin/variance controls, Stack transfer and donor-specific CD4 responses | Residual10 emission and V24 helped. V20 did not beat its V16 parent. [V16 results](V16_RESULTS.md), [V20 results](V20_RESULTS.md), [V24 results](V24_RESULTS.md). |
| V25-V29 | Intervention geometry, neural count decoder and control-preserving transport | No candidate passed its release checks. Better null realism alone did not ensure treatment accuracy. [V25](V25_RESULTS.md), [V26](V26_RESULTS.md), [V27](V27_RESULTS.md), [V28](V28_RESULTS.md), [V29](V29_RESULTS.md). |
| V30-V34 | STRING/construct features, guide correction, paired transfer, tabular foundation model and flow model | Mean/shuffled-control requirements failed. [V30](V30_RESULTS.md), [V31](V31_RESULTS.md), [V32-V34](V32_V34_RESULTS.md). |
| V35-V37 | Experimental-response retrieval, V24 residual correction, population-view training and ensemble blends | V36 released the retrieval-contrast correction; the full population-view branch was rejected. V36 and V55 half-blending lost to V24. [V35/V36](V35_V36_RESULTS.md), [V37](V37_ENSEMBLE_RESULTS.md), [V55](V55_LEADERBOARD_RESULTS.md). |
| V38-V44 | Activity/direction prediction, KOLF data, contrastive identity, LINCS fingerprints and extra residual dispersion | Failed relevant gates or target-shuffled controls. [Combined results](V38_V44_RESULTS.md). |
| V45-V53 | Residual reduction, population mixtures, corroboration, margin projection/rounding, response priors and joint mean/rank networks | Numerical gains did not establish biological transfer. Published V50 lost to V24; joint heads did not justify a release. [V45-V50](V45_V50_RESULTS.md), [V52](V52_JOINT_RESULTS.md), [V53](V53_JOINT_GUARD_RESULTS.md). |
| V56-V58 PIE | Genuine checkpoint plus saved training-response evidence; three retained outputs; 154 seen and 146 unseen validation targets; response checks before calibrated count generation | Downloaded 17.12 GiB of summaries, checked all 64,260 evidence rows, evaluated 366 eligible response cases and 176 count comparisons. Gains did not reliably beat shuffled/evidence-only controls or validate the joint decoder. Source studies used in checkpoint pretraining were not unseen-study tests. [PIE review](V56_PIE_RELEASE_REVIEW.md), [V58 results](V58_PIE_CROSS_RESULTS.md). |
| V59-V64 | State abundance/within-state response, depth adjustment, assay/common-response transformation and identity guards | State transfer failed; V61 regressed. V63 improved hosted score despite failed public gates, exposing proxy mismatch. V64 failed count checks. [State review](V59_V61_STATE_REVIEW.md), [validation audit](V62_V64_VALIDATION_AUDIT.md). |
| V65-V68 | Available-axis normalization, assay strength, null auditing and further response alternatives | V66 supplied a tiny reach-driven hosted gain; five other components declined versus V63. Null evidence limited later alternatives. [V65](V65_RESULTS.md), [V66](V66_RESULTS.md), [V67](V67_RESULTS.md), [V68](V68_NULL_AUDIT.md). |
| V69-V82 | Probability-coupled count mixture and 12 follow-up hypothesis families | V69 regressed; V70-V81 failed their original public requirements. Raw-statistic replay audited 410 completed comparisons and 3,280 tied-rank coordinates; these reuse biological cases. V82 restored identical V66 predictions. [Follow-up results](POST_V69_RESULTS.md). |
| V83-V84 | New essential-screen coverage and paired response-transfer operator; measured-mean diagnostic | New screen had zero direct target overlap. Operator trained on 638 shared interventions (513 fit/125 held out) but failed public gates. Supplying measured treatment means improved a non-deployable diagnostic markedly; these oracle labels cannot be used for final prediction. [Research](POST_V82_RESEARCH.md), [coverage](V83_COVERAGE_AUDIT.md), [V84](V84_RESULTS.md). |
| V85 | Split-sample guide-dose contrasts with permuted-guide null | Corrected sparse-log bias; no source target met both fixed dose requirements. Exact abstention; no duplicate submission. [Guide contrast](V85_GUIDE_CONTRAST.md). |
| V86 | Covariance-mixture Bayesian response prior, held-study exclusion and interval projection | Small profile gains on 655 reused pairs; shuffled control competitive and original requirements failed. Finite 100 updates did not converge. [V86](V86_RESULTS.md). |
| V87 BOLERO | Reviewed geneRegPred paper, code and checkpoint requirements; independent conditional low-rank response operator | Full DNA foundation model was not run. Compact adaptation failed constant/shuffled controls; no artifact. [Paper/code review](BOLERO_REVIEW_V87.md), [results](V87_RESULTS.md). |
| V88-V90 | Converged simplex MAP, exact mixture-CDF intervals, second-seed retest and conservative V86 projection | V88/V89 failed original gates and were not released. V90 was explicitly post hoc diagnostic selection of the unchanged interval25 method; its eventual small hosted gain does not retroactively pass failed gates. [V88](V88_RESULTS.md), [V89](V89_RESULTS.md), [V90](V90_RESULTS.md). |
| Other literature and cross-domain ideas | Arc atlas/upstream commits, cross-modal foundations, response geometry, state/distribution transport and later paired-transfer work | Reviewed with controls and feasibility distinctions. Full foundation methods were not claimed as reproduced when only compact adaptations ran. [Cross-modal review](CROSS_MODAL_MODEL_REVIEW.md), [research gaps](VCC2026_RESEARCH_AND_GAPS.md), [later literature](POST_V82_RESEARCH.md). |

The central unresolved gap is predicting context-dependent downstream treatment means and DE sets in unseen biological systems. More neural layers, broader clouds, exact moment fitting or combining scalar metric winners do not by themselves supply that information. Public sources are limited, reused and assay-mismatched; seeds and correlated CD4 conditions do not create independent biological lineages. Retrospective method correlations are descriptive rather than causal.

## 5. Final V91: work completed

### Inputs, sources and missingness

The final manifest and control counts were fully audited. Of 300 final targets, **129 overlap validation and 171 are new**. Actual D/E/F predictions were rebuilt; A/B/C arrays were not relabelled.

The final public response cache has **1,487 target/context rows: 638 reused bitwise and 849 newly extracted**. Family row counts are K562 260, CD4 625, H1 28, HCT116 282 and HEK293T 292. CD4 includes three correlated conditions. There is eligible evidence for **298/300 targets**, including 169/171 new targets. **ANG and DPY30** lack eligible source responses: downstream effects remain at the control baseline, with the fixed 0.1 on-target factor still an assumption. Missing measurements stay masked.

The pinned full Orion release supplied 30,795 HCT116 cells across 109 batches and 42,155 HEK293T cells across 223 batches, including controls. Gene mapping, source offsets and full library totals were checked against the original sampler. The H1 extension supplied 8,463 treated cells for 12 new targets across public splits; shared controls were included once. CD4 merged and disjoint-donor summaries were extended under their original pinned identities. Verified cached-page concurrency sped collection without changing source samples or values.

Descriptive finite-value Spearman basal similarities were D-CD4 Stim48hr 0.850588, E-H1 0.776744 and F-H1 0.802042. These are not cell-line identifications or response-accuracy tests.

### Fixed model and count generation

1. Pool measured public responses using the established V13/V16 shrinkage/missingness logic.
2. Apply V24 `cdf50_de` donor-supported CD4 correction.
3. Apply V66 `assay75` correction with the existing 80% V24 identity-margin guard.
4. Apply V90 `interval25`: quarter-strength normal-moment interval projection, maximum correction 0.125 LFC, with the existing 80% V66 identity guard.
5. Emit raw counts using the unchanged residual10 generator and original D/E/F control templates/library totals. The seed is 20261011; the direct on-target factor is 0.1, applied once.

Source-panel centering, rank axes, CD4 teacher and covariance priors were refit for final targets. Strengths, identity guards and emission choices were not searched against final labels or scores. Both real and shuffled priors retain the original 100-update budget; convergence and interval coverage are not established. Saved responses retain V24, V66, interval25, posterior25, shared25 and shuffled diagnostic outputs; only the fixed interval25 output was emitted.

### Completed checks and artifact

- CLI upgraded from 0.2.2 to **0.3.0**; final-manifest smoke round-trip passed.
- All nine archived V24/V66/V90 response arrays replayed bitwise under the working Python 3.14 runtime; three archived full 400-cell V90 groups regenerated exactly.
- Source/input audits checked final axes, every raw control count, missingness, own-target exclusions, coverage and original source hashes.
- Numerical likelihood parity, mixture weights, finite outputs, correction bounds and both identity guards passed.
- Twelve final 400-cell preflight groups passed before full generation; all 900 groups completed.
- Lossless memory-safe staging compared every stored count, index and pointer. Official prep passed.
- Independent validation compared the full packaged payload, all 360,000 original selected cell identities/library totals and twelve fully regenerated final groups.
- Repository audit passed for 62 frozen code files, six final inputs, source receipts and historical-protocol preservation. Frozen legacy helper line endings were preserved explicitly; no scientific values changed.

**File:** `experiments/final_v91/prediction.vcc`

**Bytes:** 3,929,180,160

**SHA256:** `7595e024a9a4ceac08937e0b673dd056be58bb005afc8fbfe2b7f35f7565bb88`

**Panel:** `test / vcc2026-test-1`; **360,000 cells**, **18,533 genes**, **2,283,563,644 nonzeros**; maximum cell total **53,865**.

These checks establish construction, provenance and reproducibility. They do not establish final treatment accuracy. The final controls contain no treatment truth; the local validations cannot certify a ranking improvement.

## 6. Reproduction and evidence map

Root code is versioned; large controls, raw sources, model arrays, H5ADs and VCC files remain local and ignored. A Git clone alone does not contain the large data; consult pinned-source collectors and audit manifests. Use `.competition_env\Scripts\python.exe -X utf8` from the root. The exact final runtime dependencies are in [requirements-final.txt](requirements-final.txt).

| Stage | Root code | Main evidence under `experiments/final_v91/analysis/` |
|---|---|---|
| Final input inspection | `inspect_final_inputs.py` | `input_audit.json` |
| Public data acquisition | `collect_final_cd4.py`, `collect_final_orion.py`, `collect_final_h1.py` | `cd4_collection.json`, `HCT116_new_source.json`, `HEK293T_new_source.json`, `H1_new_source.json` |
| Source assembly and coverage | `build_final_sources.py`, `audit_final_sources.py` | `source_provenance.json`, `source_audit.json`, `target_coverage.csv` |
| Runtime and CLI compatibility | `audit_final_replay.py`, `test_final_count_replay.py`, `test_final_cli.py` | `validation_replay_A/B/C.json`, `count_runtime_replay.json`, `cli_final_smoke.json` |
| Basal diagnostics | `diagnose_final_contexts.py` | `basal_rank_correlations.csv`, `context_diagnostic.json` |
| Freeze, fit, retain controls | `final_response_methods.py`, `fit_final_submission.py`, `audit_final_model.py` | `plan.json`, `reference_model.json`, `model_audit.json`, model metadata |
| Preflight and full counts | `build_final_submission.py`, `summarize_final_generation.py` | `prediction_validation.json`, `all_group_moment_summary.json`, group statistics |
| Official packaging and independent validation | `package_final_submission.py`, `validate_final_submission.py` | `memory_safe_stage.json`, `submission_ready.json`, `independent_validation.json` |
| Ordered orchestration/reporting | `run_final_pipeline.py`, `report_final_results.py`, `audit_final_repository.py` | `experiment_complete.json`, `repository_audit.json` |
| This dated documentation | `document_project_progress.py`, then `--check` | `analysis/project_progress_20261023/` |

The ordered pipeline intentionally never uploads or pushes. **Do not rerun generation, fitting or packaging against the active uploading file.** Detailed preparation and stage arguments are in [FINAL_ROUND.md](FINAL_ROUND.md); file/readiness and the original CMD command are in [FINAL_V91_RESULTS.md](FINAL_V91_RESULTS.md). Runtime replays establish technical compatibility only.

Local scientific/release commits completed before this documentation: `bbac302` (audited final source transfer and V90 publication), `4c29c09` (source freeze and unsupported targets), `689b3e8` (validated final release). They were committed locally; no GitHub push occurred automatically. The documentation commit follows these commits.

## 7. Upload status and remaining evidence

As of the user's **23 October 2026** message, **Final V91 is uploading**. No entry ID, server completion receipt, score or rank has been supplied. This is a user-reported state, not an authenticated server query. The [separate status snapshot](analysis/project_progress_20261023/upload_status.json) records that distinction. Frozen readiness/experiment-completion JSON retains `uploaded: false`, describing the earlier local preparation stage.

The remaining action is to let the current upload finish and retain its entry ID and completion receipt. Confirm the server receipt identifies the test partition/final panel. If interrupted, use the CLI resume flow for the existing entry instead of creating a duplicate. Preserve the same package while uploading. Final scores/rank stay blind until the organizer publishes them.

Only the team's latest final entry counts. The supplied deadline is **5 November 2026, 23:59 UTC** (**6 November, 00:59 Europe/Oslo**); the allowance is two submissions per UTC day. See the [official CLI guide](https://vcc-cli-wiki.virtualcellchallenge.org/). No extra upload, cancellation, server status query or GitHub push is performed by this documentation task.
