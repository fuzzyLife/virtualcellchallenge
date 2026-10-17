# V52 joint mean/rank training results

Shared neural predictor with separate mean log2FC and Mann-Whitney rank-effect outputs. [Frozen protocol and reproduction](V52_JOINT_STRATEGY.md). Public data are reused; D3 remains unopened. Mean/rank arrays require separate count-decoder validation before a competition release.

## Extracted paired-view labels

| Source | Study | Supported targets | Eligible target/gene pairs |
|---|---|---:|---:|
| H1 | H1 | 12 | 128,357 |
| HCT116 | Orion | 139 | 1,404,204 |
| HEK293T | Orion | pending | pending |
| CD4_Rest_D2 | CD4 | pending | pending |
| CD4_Stim48_D2 | CD4 | pending | pending |

The two views are sampling replicates, not independent studies. Missing target/gene labels abstain. All measured genes were considered; each fitting target samples at most 768 eligible gene pairs, while outer evaluation uses every eligible measured gene.

## Fitted models

| Excluded study | Model | Selected epoch | Training pairs | Source-validation pairs |
|---|---|---:|---:|---:|

Whole-study evaluation is pending.

## Integrity and delivery

- label_audit.json: pending
- boundary_audit.json: pending
- final_audit.json: pending

Model files live under experiments/v52_joint_mean_rank/analysis/models. The all-source joint.pt checkpoint is an exploratory inference artifact; it does not have an independent all-source test. Full-axis held-study joint mean/rank predictions are under analysis/predictions. Each model and prediction matrix has a SHA256 sidecar. Large binaries are ignored by Git.

Root API: predict_v52_joint_response.predict(query_cpm, outer="all") returns mean effects, rank effects and a supported-coordinate mask. It keeps V24 mean effects for output genes below the rank-label expression gate and excludes the own gene. A rank effect is not a p-value. Current source-control features are means/ranks and donor response statistics, not a fitted full query-control distribution; a later rank-aware count decoder would need that distribution and tie structure. Production inference has more available donor studies than some fitting folds; this requires further count and biological validation.

Audits: independent raw-cell pairwise rank replay; treatment/control view disjointness; excluded-response mutation; source-validation target-group boundaries; model and inference array hashes/bounds. These checks verify the implementation and exclusions, not hidden-context generalization. No new .vcc file was created or submitted.

~~~bat
SCIENCE_PYTHON audit_v52_joint_response.py labels
SCIENCE_PYTHON audit_v52_joint_response.py boundaries
SCIENCE_PYTHON audit_v52_joint_response.py artifacts
SCIENCE_PYTHON document_v52_joint_response.py
~~~
