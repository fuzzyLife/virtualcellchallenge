# Joint mean/rank predictor: V52 and corrected V53

Implemented and trained a shared 128/64 neural network for mean log2FC and rank effect relative to controls, with three whole-study outer exclusions, single-head ablations, linear heads and a trained shuffled-target control. All-source checkpoints use the same fixed fitting schedule. [V52 protocol](V52_JOINT_STRATEGY.md), [feature-contract correction](V53_JOINT_GUARD_STRATEGY.md), [original V52 results](V52_JOINT_RESULTS.md).

V53 reuses exact V52 labels and changes only how fitting-constant features are handled. It was frozen after partial predecessor evidence; it is not a blind new hypothesis. Both versions use repeatedly examined public biology; D3 remains unopened.

## Equal-study held-out errors

| Version | Predictor | Mean MSE | Rank MSE |
|---|---|---:|---:|
| V52 | v24_pooled_rank | 0.1088581 | 0.0016190 |
| V52 | zero_rank | 0.1088581 | 0.0011819 |
| V52 | linear | 0.1085070 | 0.0019739 |
| V52 | joint | 0.1198082 | 0.0016215 |
| V52 | mean_only | 0.1265976 | 0.0016190 |
| V52 | rank_only | 0.1088581 | 0.0018185 |
| V52 | shuffled_joint | 0.1251375 | 0.0015417 |
| V53 guarded | v24_pooled_rank | 0.1088581 | 0.0016190 |
| V53 guarded | zero_rank | 0.1088581 | 0.0011819 |
| V53 guarded | linear | 0.1085041 | 0.0019739 |
| V53 guarded | joint | 0.1106324 | 0.0015644 |
| V53 guarded | mean_only | 0.1109860 | 0.0016190 |
| V53 guarded | rank_only | 0.1088581 | 0.0018042 |
| V53 guarded | shuffled_joint | 0.1204421 | 0.0016098 |

Corrected joint mean MSE versus V24: **+1.63%**. Rank MSE versus pooled rank: **-3.38%**; versus zero rank: **+32.36%**. Negative error changes are favorable.

The zero-rank baseline retains V24 means and predicts theta=0; it is not a zero-expression submission. Mean/rank errors are raw response-label diagnostics, not official VCC scores. Directional and significance metrics require count construction and an independent reference.

## Prespecified pilot checks

- mean_not_worse: **False**
- rank_gain_5pct: **False**
- beats_zero_rank: **False**
- beats_shuffled_both: **True**

Pilot success: **False**.

The corrected joint model does not meet the fixed response-level checks. Keep V24 as the preferred submission; the trained checkpoint is an exploratory artifact, not an established improvement.

## Study-specific errors

| Held-out study | Joint mean / V24 MSE | Joint rank / pooled MSE | Joint rank / zero MSE |
|---|---:|---:|---:|
| H1 | 1.0309 | 0.7374 | 1.0920 |
| Orion | 1.0521 | 1.2760 | 1.3541 |
| CD4 | 0.9609 | 0.9196 | 1.6850 |

## Label quality and limitations

There are 366 supported source-target groups: H1 12, HCT116 139, HEK293T 147, CD4 Rest D2 34 and CD4 Stim48 D2 34. All measured output genes are considered; fitting uses at most 768 eligible pairs per group and held-study evaluation uses all eligible measured outputs. Sample replicates are not new biological contexts. Hash0 targets are excluded from fitting labels, but allowed cross-study experimental responses to them remain retrieval features.

The paired-view rank correlations are approximately 0.691 in H1, 0.022 in HCT116, 0.249 in HEK293T and 0.62 in CD4. This weak Orion reproducibility makes noise shrinkage an especially important negative control. The complete source/view diagnostic is in V52 analysis/label_noise.json.

The feature guard prevents unsupported train-constant donor channels from activating untrained weights. It does not establish context-specific biological transfer or eliminate all distribution shift. Source mean labels and theta are lane-matched; query features use pooled control profiles. The model does not consume a fitted full control CDF/dispersion representation. A useful future extension would need source-validated control-distribution features and evidence that it beats zero/shuffled rank predictions, before a count decoder is trusted.

## Deliverables

- Corrected all-source checkpoint: experiments/v53_joint_guard/analysis/models/all/joint.pt, with SHA256 and fitting history in joint.json.
- Corrected full-axis held-study mean/rank arrays: experiments/v53_joint_guard/analysis/predictions/{source}.npz, with source, masks and hashes in the JSON sidecars.
- Original V52 models and arrays remain under experiments/v52_joint_mean_rank/analysis.
- Root training, inference, tests and auditing scripts are versioned. Binary labels/models/arrays remain local and regenerable.

~~~python
from run_v53_joint_guard import install
install()
from predict_v52_joint_response import predict
result = predict(query_cpm, outer="all")
mean_effects = result["mean"]
rank_effects = result["rank"]
eligible = result["valid"]
~~~

query_cpm must contain 18,533 per-cell-mean CPM values in the competition gene order. The function returns all 300 targets by default; genes below the query control expression gate keep the V24 mean effect, with rank eligibility marked absent; the mask is not a guarantee of direct donor evidence. Own-target coordinates are excluded. The output is not an upload-ready count matrix.

## Audit status

- label_audit.json: **passed**
- boundary_audit.json: **passed**
- feature_contract_audit.json: **passed**
- final_audit.json: **passed**
- delivery_manifest.json: **passed**

The [delivery manifest](experiments/v53_joint_guard/analysis/delivery_manifest.json) binds the all-source checkpoint and API code to completed audits and a CPU reload/forward smoke test. No new .vcc file was created or uploaded. No Git push was performed.
