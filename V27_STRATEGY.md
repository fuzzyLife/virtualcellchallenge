# V27: query-control anchors and relative neural responses

This is a prospective experiment, not a claimed leaderboard improvement. V24 remains the published reference (overall 0.14299785521745675).

## Why these tests

V26 reduced H1 null differential-expression calls from 821 to 352 but failed calibration. Of its calls, 206 were in 1,061 fallback coordinates. Its empirical comparator produced 6,764 calls after measurement resampling and moment repairs. These findings motivate separately testing sampler stages, fallback locking, query-control quantile calibration, relative neural response changes, and learning panel-gene output variation. They do not establish that any proposed correction works.

## Frozen design

`prepare_v27_experiment.py` hashes dependencies, historical input evidence, prospective populations, methods, checks, and the sole eligible finalist before fitting or scoring. The matched network uses the V26 architecture and source sampling, 12 epochs, and held-cell validation ELBO selection. The panel model permits shared measured target-gene outputs but masks the directly perturbed coordinate per example. Inputs and state/response bases remain nonpanel. Source group means use fitting cells only. Every evaluation family is excluded from source fits and mean prediction. H1's masked checkpoint is the immutable V26 conditional fit.

Nine paired arms: V24, V26, locked fallback, raw anchor, calibrated anchor, relative neural change, calibrated relative change, calibrated panel-relative change, and target-descriptor shuffle. Only `panel_relative` is eligible; diagnostic arms cannot be promoted after seeing scores. Both learned relative arms use D(z,state,effect) minus D(z,state,zero), with common Gaussian latent draws. They do not identify paired counterfactual cells.

Query-control calibration maps normalized template order statistics to a disjoint NTC fitting population, with independent tie breaking. Templates, fitting controls, actual controls and reference controls have disjoint cell identities. Splits are cell-disjoint, not claimed guide- or batch-disjoint; overlap is explicitly archived. H1 has only 1,536 available controls, leaving 336 fitting controls for each null split. Three null splits reuse biology and are robustness probes, not independent-study validation.

Every final population retains exactly the V24 generated gene totals and template libraries. Mean CPM-probability L1 drift must be <=0.0005. Integer projections can distort distribution and are evaluated after all repairs. The isolated fallback arm locks those coordinates to V24 counts and restricts integer repairs; it is diagnostic and may fail the mean envelope. It does not test preserving raw fallback shape independently of V24.

## Execution and stopping

First replay the V26 empirical null through successive sampler stages. Fit the matched H1 panel-output model; score all three null splits and all eight H1 treatment cases across all nine arms. The treatment pilot is explicitly diagnostic even after null failure so all proposed ideas receive actual intervention tests. H1 and the public development cases are spent biology, not a new reserve. Passing zero-response identity alone cannot justify release.

Expand to HCT116, HEK293T and CD4 only if the fixed finalist passes all H1 null and treatment checks. The remaining complete design has 44 historical treatment cases. Only the locked finalist can open a new donor reserve after development passes; all CD4 evidence must be excluded for that reserve. Its partially downloaded raw source must be completed and author-MD5 verified first. No automatic upload.

Null checks retain V26 envelopes: calls <=max(20, real+10), detection MAE <=max(.03,2*real+.005), log-variance MAE <=max(.5,1.5*real+.1). Treatment checks retain historical raw guards, at least three primary improvements including DE, worst-family fidelity/reach, improved centered shape Wasserstein, and at least two DE improvements over each anchor, calibrated-anchor, masked-calibrated-relative and shuffled-descriptor comparator. HCT116 lacks NMAE and CD4 D1 has 200 actual treatment cells; these limitations cannot be disguised as complete host-metric validation.

## Reproduce

Use the Python 3.12 scientific runtime and prior verified source extractions, V24 teacher and V25 response bases. Official packaging uses a separate runtime and must not regenerate scientific arrays.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
set "VCC_SCI_PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
"%VCC_SCI_PY%" test_v27_controls.py
"%VCC_SCI_PY%" run_v27_pipeline.py
```

Frozen evidence refuses overwrite; identical verified caches may resume. Root scripts, small hash/score records and documentation are versioned. Raw count arrays, weights, caches and submission packages are regenerable and ignored by Git. No package is upload-ready without official prep, complete payload/axis/provenance checks and exact replay audits.

## Operational execution record

The frozen original pipeline remains intact. After its complete null comparisons, an operational handoff preserved all caches and resumed treatment comparisons with `run_v27_with_processes.py pipeline`. This adapter runs four independent case processes with two BLAS threads, unchanged seeds/methods/guards, and verified cached identities. `test_v27_processes.py` checks exact serial/process Atlas arithmetic; full-axis independent null/treatment replays check scientific output parity. Records are `analysis/orchestration_handoff.json` and `process_execution.json`. The additional root `audit_v27_cell_labels.py` verifies actual target/control cell labels without modifying the originally executed training-input auditor.

### Existing records and fresh preparation

The recorded checkout already contains its frozen plan. Skip the preparation command when replaying it; use the pipeline, auditor and reporter with the verified prerequisite data/arrays. `prepare_v27_experiment.py` is a once-only step for an empty new analysis directory, not an overwrite command. A fresh retraining or new workspace needs its own run directory/manifest and new fitting-time records; retain these historical receipts and validate regenerated numerical inputs against their declared hashes. Large raw/weight/count arrays are ignored by Git and must be reconstructed from the preceding source/pretraining steps or retained verified caches. Exact historical replay requires the saved weights and pinned numerical runtime.
