# V5 published leaderboard result

Entry `BPI0n52HTsG0EYTFOf9H` is **published**, with overall **-0.113637**. V4 remains the best project submission at **-0.050674**. All four published entries use the same validation panel and anchor set; V3 has no recorded score.

| Scaled metric | First | V2 | V4 | V5 | V5 minus V4 | Contribution to overall change |
|---|---:|---:|---:|---:|---:|---:|
| overall | -0.276290 | -0.223325 | -0.050674 | -0.113637 | -0.062963 |  |
| pds | 0.019644 | 0.103042 | 0.255302 | 0.447856 | +0.192554 | +0.032092 |
| mse | 0.000000 | 0.000000 | 0.000000 | 0.035294 | +0.035294 | +0.005882 |
| nmae | 0.016006 | 0.032352 | 0.067539 | 0.070189 | +0.002650 | +0.000442 |
| fid | -1.665331 | -1.470413 | -0.683855 | -1.329100 | -0.645246 | -0.107541 |
| reach | 0.033874 | 0.041380 | 0.086849 | 0.121061 | +0.034213 | +0.005702 |
| jac | -0.061932 | -0.046314 | -0.029878 | -0.027120 | +0.002758 | +0.000460 |

Five of six scaled components improve; direction fidelity changes by **-0.645246**. Its contribution is **-0.107541**, while the other five contribute **+0.044578**, leaving an overall change of **-0.062963**. This is an arithmetic decomposition of the six-component mean, not causal attribution.

| Raw metric | V4 | V5 | V5 minus V4 | Improvement |
|---|---:|---:|---:|---|
| pds_cosine | 0.616462 | 0.702367 | +0.085905 | yes |
| expr_mse_unbiased_capped_norm | 1.247733 | 0.970079 | -0.277654 | yes |
| de_wilcoxon_lfc_nmae | 0.959164 | 0.957501 | -0.001663 | yes |
| de_wilcoxon_direction_fidelity_yield_raw | 0.309235 | 0.118548 | -0.190687 | no |
| de_wilcoxon_direction_reach_raw | 0.156391 | 0.186364 | +0.029973 | yes |
| de_wilcoxon_sig_jaccard | 0.019731 | 0.020722 | +0.000991 | yes |

Raw expression error falls **22.25%**. Direction-fidelity yield falls **61.66%**. Fidelity yield is response-budget-normalized; it is not simple sign precision. The supplied aggregate metrics cannot distinguish context-specific failure, response undercoverage or incorrect directions.

V5 simultaneously changes public H1 coverage, B/C response calibration, moment matching, control templates and A promoter corrections. No component or context can be blamed from the aggregate result. A useful follow-up is controlled source-domain ablation of these changes against V4, measuring mean/pooled effects, variance, response coverage and direction yield together. The new 300-gene neural models reconstruct control cells; their reconstruction gains do not establish better CRISPRi response predictions.

Rank was **1100** in the user report; the recorded API check returned **1099**. Ranks depend on the current leaderboard and are not directly comparable across publication dates. Submitted at `2026-10-06T06:26:47.603765Z`; reported on 2026-10-07. Partition `val`, panel `vcc2026-val-1`, anchors `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`.

## Records and reproduction

The original `submission_ready.json` remains the generation-time snapshot. Publication is recorded separately in `leaderboard_result.json`, `published_api_status.json`, `leaderboard_comparison.csv` and `raw_metric_comparison.csv` under `experiments/v5_atlas_calibrated/analysis/`. The local package hash is inherited from the verified readiness record; the status API does not independently expose or confirm that hash.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc status BPI0n52HTsG0EYTFOf9H
.\.competition_env\Scripts\python.exe record_v5_result.py
```

Both commands read the existing entry; they do not upload a duplicate. [V5_SUBMISSION.md](V5_SUBMISSION.md) and [V5_MODEL_STRATEGY.md](V5_MODEL_STRATEGY.md) preserve generation and verification details.
