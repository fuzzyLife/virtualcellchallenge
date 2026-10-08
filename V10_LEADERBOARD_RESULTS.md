# V10 published leaderboard result

Published as `CYpAqJP7bUYqjA2ouD7U`, overall **0.116754121**, rank **585** in the user report and **585** at the recorded API check. V10 is the best published project submission, improving on V4 by **+0.167428032**. Ranks change over time. V4 and V10 use the same partition, panel and anchors.

## Scaled score comparison

| metric | v4 | v9 | v10 | change_from_v4 | overall_contribution_change_from_v4 |
| --- | --- | --- | --- | --- | --- |
| overall | -0.050674 | -0.070447 | 0.116754 | 0.167428 | undefined |
| pds | 0.255302 | 0.313301 | 0.572464 | 0.317162 | 0.052860 |
| mse | 0.000000 | 0.033719 | 0.049815 | 0.049815 | 0.008302 |
| nmae | 0.067539 | 0.074927 | 0.049886 | -0.017653 | -0.002942 |
| fid | -0.683855 | -0.921183 | -0.022335 | 0.661520 | 0.110253 |
| reach | 0.086849 | 0.102034 | 0.046547 | -0.040302 | -0.006717 |
| jac | -0.029878 | -0.025482 | 0.004148 | 0.034026 | 0.005671 |

Four scaled components improve versus V4: discrimination, expression accuracy, direction fidelity and significance overlap. Fold-change accuracy and reach decline. Discrimination contributes **+0.052860347** and fidelity **+0.110253279** to the overall gain. Together they account for **97.42%** of the net improvement. Contributions are arithmetic changes in the six-metric mean, not causal attribution or significance estimates. Scaled MSE/NMAE are accuracy scores: higher is better; V4's zero expression score was clipping, not zero error.

## Raw metric comparison

| metric | v4 | v10 | change_from_v4 | improved_from_v4 |
| --- | --- | --- | --- | --- |
| pds_cosine | 0.616462 | 0.759535 | 0.143073 | True |
| expr_mse_unbiased_capped_norm | 1.247733 | 0.945955 | -0.301778 | True |
| de_wilcoxon_lfc_nmae | 0.959164 | 0.970311 | 0.011147 | False |
| de_wilcoxon_direction_fidelity_yield_raw | 0.309235 | 0.505983 | 0.196748 | True |
| de_wilcoxon_direction_reach_raw | 0.156391 | 0.120318 | -0.036073 | False |
| de_wilcoxon_sig_jaccard | 0.019731 | 0.031815 | 0.012084 | True |

Raw fidelity yield rises **63.62%**, expression error falls **24.19%**, and reach falls **23.07%**. Direction fidelity remains below the anchor level (its scaled component is still negative). Raw metrics and anchored scores are different quantities.

## What the result establishes

V10 beats V4 on the scored validation panel. It changes both the response estimator (five-source pooling, bounded knockdown-dose correction and uncertainty shrinkage) and count emission (compact synthetic populations with two-moment integer fitting). The result does not isolate either change's contribution. Corrected public validation anticipated the average fidelity gain and the reach decline, but source-family behavior was uneven; HCT116 fidelity regressed with full references. Compact emission still has substantial neutral-control DE bias. A leaderboard improvement does not establish faithful biological single-cell variance.

Keep the exact V10 artifact and frozen parameters as the current published baseline. Next compare V4 responses with the same compact emitter to isolate response pooling from emission; reserve any amplitude/dispersion or coverage choices for a new experiment. Aim to recover reach and fold-change accuracy while protecting fidelity and discrimination. Do not select another model using only expression error, or infer tissue-specific causes from aggregate scores. One submission cannot quantify random-seed uncertainty or establish final-round generalization.

## Provenance and reproduction

Entry submitted `2026-10-08T17:09:31.389918Z`; user-reported upload duration 17 minutes 5 seconds. Partition `val`, panel `vcc2026-val-1`, anchors `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`. The independently verified local package SHA-256 is `0aea2ed2900a404de5b99963a9e00887a1b1fff1c14165e631869d556374da44`. The status API does not expose the uploaded package hash; local integrity and publication are recorded separately.

Generation-time upload=false snapshots, model code, choices and validation records remain unchanged. `publication_manifest.json` links those snapshots with the saved API response, result, comparisons and publication documents. Large prediction files remain excluded from Git.

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc status CYpAqJP7bUYqjA2ouD7U
.\.competition_env\Scripts\python.exe record_v10_result.py
```

For offline reproduction from the saved status response, add `--api-file experiments/v10_population_transport/analysis/published_api_status.json`. These commands inspect/document the existing entry; they do not submit.
