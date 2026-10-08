# V9 published leaderboard result

Entry `QBpef7kQNiFa4IFSxNi5` is **published**, overall **-0.070447281**, rank **1084** in the user report and **1086** at the read-only API check. Ranks change over time. V4 remains the best published project submission at **-0.050673910**. Both entries use the same partition, panel and anchors.

## Scaled score comparison

| metric | v4 | v7 | v9 | change_from_v4 | overall_contribution_change_from_v4 |
| --- | --- | --- | --- | --- | --- |
| overall | -0.050674 | -0.051387 | -0.070447 | -0.019773 | undefined |
| pds | 0.255302 | 0.235853 | 0.313301 | 0.057999 | 0.009667 |
| mse | 0.000000 | 0.000000 | 0.033719 | 0.033719 | 0.005620 |
| nmae | 0.067539 | 0.071861 | 0.074927 | 0.007388 | 0.001231 |
| fid | -0.683855 | -0.672676 | -0.921183 | -0.237328 | -0.039555 |
| reach | 0.086849 | 0.086715 | 0.102034 | 0.015186 | 0.002531 |
| jac | -0.029878 | -0.030073 | -0.025482 | 0.004396 | 0.000733 |

Five of six scaled metrics improve versus V4. Fidelity contributes **-0.039554644** to the overall difference, while the other five combined contribute **+0.019781274**. Their sum is **-0.019773371**. These are arithmetic contributions to the six-metric mean, not causal attribution or significance estimates.

## Raw metrics

| metric | v4 | v9 | change_from_v4 | improved_from_v4 |
| --- | --- | --- | --- | --- |
| pds_cosine | 0.616462 | 0.642059 | 0.025596 | True |
| expr_mse_unbiased_capped_norm | 1.247733 | 1.011257 | -0.236476 | True |
| de_wilcoxon_lfc_nmae | 0.959164 | 0.955218 | -0.003946 | True |
| de_wilcoxon_direction_fidelity_yield_raw | 0.309235 | 0.238046 | -0.071189 | False |
| de_wilcoxon_direction_reach_raw | 0.156391 | 0.169710 | 0.013319 | True |
| de_wilcoxon_sig_jaccard | 0.019731 | 0.021359 | 0.001628 | True |

Raw expression error falls **18.95%**; raw fidelity yield falls **23.02%**. The positive scaled expression score indicates improved accuracy relative to the anchors; V4's zero scaled expression score was clipping, not zero error. Raw discrimination, LFC error, reach and overlap also improve. A larger negative scaled fidelity component defeats those gains.

## What we learned

V9 preserved the exact V4/V7 biological effect arrays and added 50% mean/pooled-profile correction on top of V4 zero-aware simulation, using new random draws. This is a decoder experiment, not a newly learned biological response model. The result demonstrates that improved mean expression, discrimination and neutral-control preservation do not establish perturbation direction fidelity. Aggregate scores cannot show whether the loss comes from fewer correct calls, extra incorrect calls, changed ranks/variance, or specific contexts/targets. A single upload cannot estimate seed uncertainty or isolate each emission change causally.

The public confirmation proxy improved by +0.019533 and passed its declared guard. However, public mean raw fidelity already fell from **0.479674** to **0.397573**. The guard used a reference-call>=10 budget-weighted subset, whose fidelity changed from 0.002476 to 0.002724. That very low, narrowly supported improvement did not protect the overall yield. Missing metric panels also reduced the effective support of public utility averages. The recorded confirmation acceptance remains historically true under that declared rule, but it did not establish hidden improvement.

A retrospective full mean-fidelity guard allowing at most 0.01 loss would have **rejected** V9 on public confirmation. This is an audit finding, not a preregistered rule or a retuned confirmation success. The original selection and confirmation records remain unchanged.

## Next decision

Keep V4 as the published reference. Retain corrected common-control scoring and expression jackknifes, but do not promote V9's count decoder. Before another upload, require fidelity protection both on a common supported-target set and on overall yield with explicit missingness; report real/predicted/correct DE budgets per target and equal-lineage/context results. A gain confined to the weighted subset should not override a material overall fidelity loss. Add neutral and perturbation distribution/detection checks together, and use paired seeds with reference DE budgets closer to the challenge. New thresholds/folds belong to a new experiment; V9 confirmation is now used development data.

Prioritize validated response direction and target-specific uncertainty over another unchanged-model resampling upload. Do not claim a particular context, dose or gene is responsible from these aggregate numbers. This result alone does not justify uploading a 25% decoder variant or predict that it will score better.

## Provenance and reproduce

Submitted `2026-10-08T11:51:33.789090Z`; upload duration 13 minutes 20 seconds from the user report. Partition `val`, panel `vcc2026-val-1`, anchors `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`.

The independently validated local package SHA-256 is `b16504d89bf238d9b6fa3c38ba65a8108adb43eab11cbd9c556e853ecf15f8cf`. The API does not expose the uploaded package hash; local identity and publication are recorded separately. Generation-time upload=false snapshots, code/data hashes and calibration decisions are preserved. `publication_manifest.json` records their hashes and links the publication overlay, status JSON and scaled/raw comparison CSVs.

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc status QBpef7kQNiFa4IFSxNi5
.\.competition_env\Scripts\python.exe record_v9_result.py
```

For an offline reproduction from the saved status JSON, add `--api-file experiments/v9_activation_preserving/analysis/published_api_status.json`. These commands inspect/document the existing entry; they do not upload. [V9_SUBMISSION.md](V9_SUBMISSION.md) retains the original numerical strategy and public evaluation results.
