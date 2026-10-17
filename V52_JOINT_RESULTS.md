# V52 joint mean/rank training results

Shared neural predictor with separate mean log2FC and Mann-Whitney rank-effect outputs. [Frozen protocol and reproduction](V52_JOINT_STRATEGY.md). Public data are reused; D3 remains unopened. Mean/rank arrays require separate count-decoder validation before a competition release.

## Extracted paired-view labels

| Source | Study | Supported targets | Eligible target/gene pairs |
|---|---|---:|---:|
| H1 | H1 | 12 | 128,357 |
| HCT116 | Orion | 139 | 1,404,204 |
| HEK293T | Orion | 147 | 1,594,407 |
| CD4_Rest_D2 | CD4 | 34 | 311,770 |
| CD4_Stim48_D2 | CD4 | 34 | 319,554 |

### Paired-view label reproducibility

| Source | Mean correlation | Rank correlation | Target-centered rank correlation |
|---|---:|---:|---:|
| H1 | 0.666 | 0.691 | 0.653 |
| HCT116 | 0.055 | 0.022 | 0.013 |
| HEK293T | 0.118 | 0.249 | 0.109 |
| CD4_Rest_D2 | 0.538 | 0.622 | 0.613 |
| CD4_Stim48_D2 | 0.528 | 0.616 | 0.608 |

Weak HCT116 paired-view reproducibility means lower prediction error can reflect shrinkage of sampling noise. Zero-rank and target-shuffled controls are required to distinguish this from response learning. Centering is a separate descriptive diagnostic, never a fitted feature; reused controls across targets mean these are not independent confidence bounds.

The two views are sampling replicates, not independent studies. Missing target/gene labels abstain. All measured genes were considered; each fitting target samples at most 768 eligible gene pairs, while outer evaluation uses every eligible measured gene.

## Fitted models

| Excluded study | Model | Selected epoch | Training pairs | Source-validation pairs |
|---|---|---:|---:|---:|
| H1 | joint | 4 | 221,952 | 49,920 |
| H1 | mean_only | 5 | 221,952 | 49,920 |
| H1 | rank_only | 4 | 221,952 | 49,920 |
| H1 | shuffled_joint | 8 | 219,023 | 49,920 |
| Orion | joint | 1 | 46,848 | 14,592 |
| Orion | mean_only | 3 | 46,848 | 14,592 |
| Orion | rank_only | 1 | 46,848 | 14,592 |
| Orion | shuffled_joint | 1 | 46,646 | 14,592 |
| CD4 | joint | 7 | 190,464 | 38,400 |
| CD4 | mean_only | 8 | 190,464 | 38,400 |
| CD4 | rank_only | 7 | 190,464 | 38,400 |
| CD4 | shuffled_joint | 8 | 187,691 | 38,400 |
| all | joint | 6 | 229,632 | 51,456 |

## Whole-study raw prediction errors

| Predictor | Mean MSE | Rank MSE | Mean MAE | Rank MAE |
|---|---:|---:|---:|---:|
| v24_pooled_rank | 0.108858 | 0.001619 | 0.195921 | 0.028605 |
| zero_rank | 0.108858 | 0.001182 | 0.195921 | 0.023528 |
| linear | 0.108507 | 0.001974 | 0.210460 | 0.033039 |
| joint | 0.119808 | 0.001621 | 0.223357 | 0.029357 |
| mean_only | 0.126598 | 0.001619 | 0.231454 | 0.028605 |
| rank_only | 0.108858 | 0.001818 | 0.195921 | 0.031369 |
| shuffled_joint | 0.125137 | 0.001542 | 0.231504 | 0.028358 |

Prespecified pilot checks:
- mean_not_worse: **False**
- rank_gain_5pct: **False**
- beats_zero_rank: **False**
- beats_shuffled_both: **False**

Pilot success: **False**. This is not a private score or an upload decision.

## Per-study joint versus baseline

| Study | Joint mean MSE / V24 | Joint rank MSE / pooled | Joint rank MSE / zero |
|---|---:|---:|---:|
| H1 | 1.4901 | 0.8273 | 1.2252 |
| Orion | 1.0734 | 1.2622 | 1.3395 |
| CD4 | 1.0247 | 0.9348 | 1.7128 |

Joint changes versus the reference: mean MSE +10.06%, rank MSE +0.15%. Negative error changes are favorable. View disagreement, few independent studies, assay differences and the bounded residual architecture limit the interpretation. The same raw prediction loss does not ensure better significance, FID, REACH or a private rank.

## Integrity and delivery

- label_audit.json: **passed**
- boundary_audit.json: **passed**
- final_audit.json: **passed**

Model files live under experiments/v52_joint_mean_rank/analysis/models. The all-source joint.pt checkpoint is an exploratory inference artifact; it does not have an independent all-source test. Full-axis held-study joint mean/rank predictions are under analysis/predictions. Each model and prediction matrix has a SHA256 sidecar. Large binaries are ignored by Git.

Root API: predict_v52_joint_response.predict(query_cpm, outer="all") returns mean effects, rank effects and a supported-coordinate mask. It keeps V24 mean effects for output genes below the rank-label expression gate and excludes the own gene. A rank effect is not a p-value. Current source-control features are means/ranks and donor response statistics, not a fitted full query-control distribution; a later rank-aware count decoder would need that distribution and tie structure. Production inference has more available donor studies than some fitting folds; this requires further count and biological validation.

Audits: independent raw-cell pairwise rank replay; treatment/control view disjointness; excluded-response mutation; source-validation target-group boundaries; model and inference array hashes/bounds. These checks verify the implementation and exclusions, not hidden-context generalization. No new .vcc file was created or submitted.

~~~bat
SCIENCE_PYTHON audit_v52_joint_response.py labels
SCIENCE_PYTHON audit_v52_joint_response.py boundaries
SCIENCE_PYTHON audit_v52_joint_response.py artifacts
SCIENCE_PYTHON document_v52_joint_response.py
~~~
