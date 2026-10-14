# V25 measured results: intervention-aligned geometry

**Retain V24.** All three predeclared finalists fail development. No production package was generated, no upload was attempted, and the new-donor reserve was not scored. This is a negative result for these bounded predictive models, not a proof that every geometry-based model fails.

The frozen plan SHA256 is `883d51adb203f90e30a2473b728e986ba5961a1668394eb81efbe47da9b5943b`. Four held-family fits and one all-source fit use64 training-intervention directions on4096 shared output genes. Development emits396 populations:44 public cases times9 arms, including eight200-cell CD4RestD1 cases and36 actual400-cell cases. Baseline is the V24 algorithm under the same source-family exclusions, not the published ABC package. The biology is reused development; two CD4 conditions share a study/family.

## Family-balanced raw metrics

These are local raw metrics, not the hosted anchor-scaled leaderboard components. Higher is better for the first four; lower NMAE is better. HCT116 NMAE is unavailable, so the NMAE column averages available families only.

|Arm|PDS|Fidelity|Reach|Jaccard|NMAE|
|---|---:|---:|---:|---:|---:|
|v24|0.734524|0.572604|0.303671|0.062501|0.910671|
|force25|0.724524|0.574611|0.290197|0.062425|0.915739|
|force50|0.727857|0.570253|0.323943|0.062102|0.915300|
|bures25|0.734524|0.571648|0.267373|0.062518|0.913136|
|diagonal25|0.721190|0.573108|0.302897|0.062483|0.913793|
|shuffled25|0.731190|0.572608|0.276648|0.062304|0.912374|
|shuffled50|0.731190|0.570328|0.290373|0.062089|0.910761|
|amplitude25|0.734524|0.574139|0.267798|0.062606|0.914331|
|source_program25|0.715952|0.573458|0.287419|0.062377|0.896311|

Expression comparisons use paired per-panel ratios because simple pooled ratios are dominated by small control denominators. The mean raw expression ratio in development_summary.csv is not a hosted score. The primary expression guard passes, but none of the finalists improves three available primary metrics. Shape Wasserstein diagnostics are exploratory; this model changes response means with an unchanged residual10 emitter.

## Frozen release failures

- **force25**: three_metrics_improve, family_reach.
- **force50**: three_metrics_improve, beats_amplitude25, beats_source_program25.
- **bures25**: reach, three_metrics_improve, family_reach, beats_alpha_matched_target_shuffle, beats_amplitude25, beats_source_program25.

Force25 gains0.002008 fidelity but loses0.013475 reach and0.010000 PDS; its worst-family reach guard fails. Force50 gains0.020272 reach, while fidelity drops0.002351, PDS drops0.006667 and NMAE increases0.004629. Bures25 loses0.036298 reach and fails its shuffle/amplitude/program comparisons. The source-program control improves NMAE but loses discrimination and reach; it is a diagnostic arm, not an eligible alternative selected after seeing outcomes.

All original thresholds were retained. The conservative release_guard_amendment.json was written before opening outcome tables and can only reject further candidates. Because no finalist passes the original guards, its extra alpha0.5/Bures-specific diagnostics were not triggered. No fallback, post-hoc blend, threshold relaxation or re-selection was performed. A future experiment must freeze a new plan and retain this negative result.

## Independent verification

- Three mathematical/mechanism unit tests pass, including equal stationary covariance/Fisher geometry with different intervention responses.
- All frozen direct code/source hashes and the transitive V24 assessor hash pass.
- All five model axes, family exclusions,64-dimensional orthonormal bases, finite source shifts/weights and bounded corrections pass.
- 20 source covariances are independently recomputed from dense logCP10k counts with within-batch centering and log-depth removal.
- All 44 case identities, actual/template uniqueness, disjoint NTC halves, cache hashes and saved count digests pass.
- 15 full18533-gene count replays match exactly, with exact template libraries; 5 no-effect checks match the unchanged V16 emitter. This is not literal real-control identity or null-DE calibration.
- Recomputed release comparisons match the saved decision exactly.

The first auditor attempt incorrectly equated projected HCT116 row sums with full-source library sizes. It stopped; the auditor was corrected to use verified full-source depths, with projected sums bounded above by those depths. No frozen predictor, fitted array, test population or threshold changed. The final audit passes. Source summaries and raw CD4 covariance controls are matched at condition/study-family level, not proven exact donor aliases; whole-CD4 exclusion prevents those related sources entering CD4 validation.

## New source evidence and reserve status

Author-release metadata finds79 fibroblast and34 RPE1 competition-target overlaps, with no overlapping target reaching200 or400 retained cells. Literal symbols were checked, so low overlap/power is not a namespace conversion issue. These are CRISPRa screens and were not treated as CRISPRi labels. They remain possible representation sources requiring separate cross-modality validation. Relevant primary sources and the poster-claim assessment are in V25_STRATEGY.md.

Marson D3Rest metadata has37 overlapping targets,14 with>=400 cells. The new raw file download was stopped after development rejected every finalist: 945,815,552 partial bytes in69 range files remain resumable. The full file has **not** been MD5-verified, extracted or independently row-audited; its reserve is unopened. Extraction/audit scripts are prepared, not claimed executed. No query treatment-expression metrics were used to choose this source.

## What this changes in our reasoning

Control co-fluctuation structure is insufficient as the only context-transfer rule, even when response directions come from interventions. A bigger correction trades reach against direction/amplitude error. A source-program projection can reduce fold-change error without improving discrimination, so it cannot be promoted on one favorable metric. Future tests should identify which source intervention programs transfer with context and assay, and separate mean-response error from generated-cell DE calibration. The supplied universal-isomorphism claim supplies no verified shortcut around those measurements.

## Reproduce and inspect

Use the pinned Python3.12 scientific runtime and requirements-v23-science.txt with two BLAS threads. V25_STRATEGY.md gives from-scratch commands. Existing plans, models and decision files intentionally refuse overwrite. For the saved run:

```cmd
cd /d L:\promec\Animesh\Download\Virtual
set "VCC_SCI_PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
"%VCC_SCI_PY%" test_v25_geometry.py
"%VCC_SCI_PY%" audit_v25_experiment.py
"%VCC_SCI_PY%" document_v25_results.py
```

Numeric tables and hash evidence are in experiments/v25_cross_domain/analysis. No V25 host score is claimed. Current published champion:81Y2yX7utdAau9rr3jCj, overall0.14299785521745675.
