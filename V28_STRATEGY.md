# V28 strategy: global transport without depth-stratum repair

This separately frozen follow-up preserves V27's negative control results and ongoing matched treatment comparisons. Its sole eligible candidate is raw query-control shape transported to the family-excluded V24 means, with global gene-total/cell-library constraints and bounded CPM swaps. It has no extra multinomial measurement draw and no depth-stratum repair.

## Evidence motivating the hypothesis

On the spent first H1 null, global integer repair produces0 DE calls but CPM-probability L1 drift0.003188. Adding bounded CPM swaps without depth strata produces2 calls, detection MAE0.02843, log-variance MAE0.31690 and L1 drift0.000001196. Adding depth strata instead creates239 calls, and final swaps leave240. These are one-case numerical diagnostics, not independent validation or evidence of treatment accuracy.

The V28 plan was frozen before its treatment scores. V27 comparative treatment outcomes do not select the V28 finalist. Root `prepare_v28_experiment.py` stores hashes and exact populations in `experiments/v28_global_transport/analysis/plan.json` (SHA256 `5cf61e0d9b0031316fa01178f0646a59247d84b398c562d066c5030f7559fa20`). Historical plans and numeric code remain unchanged.

## Experiments

All three real-control null splits must pass in H1, HCT116, HEK293T, resting CD4 D1 and stimulated CD4 D2 before any V28 treatment score is opened. Fitting controls, templates, actual controls and reference controls are cell-disjoint within a case. Reused seeds/cells across cases are not independent biological replicates. Guide/batch overlap is archived, not asserted absent.

The H1 pilot has seven matched arms: V24, global raw anchor, calibrated anchor, masked relative neural change, calibrated masked relative change, calibrated panel relative change, and target-descriptor shuffle. Pretrained masked/panel weights are immutable and entirely exclude H1. The panel fit changes output masking, while its state/posterior/response features remain nonpanel; directly perturbed output is masked per fitting example. The external V25 response basis used allowed source summaries including source held-cell validation observations, so checkpoint ELBO is not wholly independent of preprocessing. Whole evaluation-family exclusion remains enforced.

Other contexts compare only V24 and the sole global-anchor finalist; no unnecessary neural fitting or state features are used. Every treatment mean predictor excludes its whole evaluation family, and inference uses only query NTCs. The full prospective development set has44 historical treatment cases (CD4 D1 has200 actual cells; others400). Intermediate context checks protect historical raw guard limits. Final context-balanced assessment additionally requires at least three primary improvements including DE, improved centered shape Wasserstein, worst-context fidelity/reach and the mean contract. HCT116 NMAE is unavailable and reported explicitly.

## Contracts and release

Every output has400 cells and18533 integer gene coordinates, identical template libraries and exact V24 generated gene totals. CPM-probability L1 drift must be <=0.0005. Global repair can still change shape; final counts, not intermediate expectations, are scored. Null thresholds are unchanged: calls <=max(20, real+10), detection MAE <=max(.03,2*real+.005), log-variance MAE <=max(.5,1.5*real+.1).

Only the raw global-anchor candidate is eligible. Neural/calibration diagnostics cannot replace it after outcomes. A failure stops further expansion and preserves V24. Only after full development passes can the locked finalist open the new D3Rest donor reserve, excluding ALL CD4 response evidence. The partial raw file must first be completed and author-MD5 verified. Passing the reserve is followed by official packaging and complete payload/axis/provenance/replay audits; no automatic upload.

This candidate transports control variation and existing predicted means. It does not claim identified target-specific stochastic dynamics. Null calibration alone cannot establish a better leaderboard score.

## Reproduce

Use the pinned scientific Python3.12 environment and prior verified source extractions/summary banks/teacher/bases. Official CLI packaging must not regenerate numerical arrays.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
set "VCC_SCI_PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
"%VCC_SCI_PY%" test_v28_transport.py
"%VCC_SCI_PY%" run_v28_with_processes.py pipeline
"%VCC_SCI_PY%" audit_v28_experiment.py
"%VCC_SCI_PY%" document_v28_results.py
```

The operational process adapter runs independent cases in four workers with two BLAS threads each. Frozen seeds, arithmetic, identities, scorer and release logic are unchanged; cached cases are verified before reuse. The shared Atlas serial/process arithmetic parity test is root `test_v27_processes.py`. Plans/results refuse overwrite; root code and small evidence are versioned, large count/model caches ignored.

### Existing records and fresh preparation

The recorded checkout already contains its frozen plan. Skip the preparation command when replaying it; use the pipeline, auditor and reporter with the verified prerequisite data/arrays. `prepare_v28_experiment.py` is a once-only step for an empty new analysis directory, not an overwrite command. A fresh retraining or new workspace needs its own run directory/manifest and new fitting-time records; retain these historical receipts and validate regenerated numerical inputs against their declared hashes. Large raw/weight/count arrays are ignored by Git and must be reconstructed from the preceding source/pretraining steps or retained verified caches. Exact historical replay requires the saved weights and pinned numerical runtime.
