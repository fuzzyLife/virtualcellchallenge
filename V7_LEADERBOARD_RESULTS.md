# V7 published leaderboard result

Entry `k0pkYuJRcw3C0U8K0NgL` is **published**, overall **-0.051386751**, rank 1052 in the user report. V4 remains best at **-0.050673910**. All six published project entries share the same validation panel and anchors.

| Scaled metric | V4 | V6 | V7 | V7 minus V4 | Overall contribution |
| --- | ---: | ---: | ---: | ---: | ---: |
| overall | -0.050673910 | -0.056790139 | -0.051386751 | -0.000712840 |  |
| pds | 0.255302018 | 0.259374252 | 0.235852667 | -0.019449351 | -0.003241558 |
| mse | 0.000000000 | 0.000000000 | 0.000000000 | +0.000000000 | +0.000000000 |
| nmae | 0.067538713 | 0.076131072 | 0.071860926 | +0.004322214 | +0.000720369 |
| fid | -0.683854648 | -0.730254749 | -0.672675846 | +0.011178802 | +0.001863134 |
| reach | 0.086848590 | 0.083530789 | 0.086714693 | -0.000133897 | -0.000022316 |
| jac | -0.029878133 | -0.029522202 | -0.030072943 | -0.000194810 | -0.000032468 |

Overall is **-0.000712840** versus V4 and **+0.005403389** versus V6. Fidelity and LFC accuracy improve; discrimination, reach and overlap decline; scaled expression remains zero. The discrimination loss offsets the fidelity and LFC gains. These contributions are arithmetic attribution within the six-component mean.

| Raw metric | V4 | V6 | V7 | V7 minus V4 | Improved versus V4 |
| --- | ---: | ---: | ---: | ---: | --- |
| pds_cosine | 0.616462282 | 0.618078781 | 0.606807878 | -0.009654404 | no |
| expr_mse_unbiased_capped_norm | 1.247732911 | 1.253227965 | 1.245516113 | -0.002216798 | yes |
| de_wilcoxon_lfc_nmae | 0.959164347 | 0.954592812 | 0.957271242 | -0.001893105 | yes |
| de_wilcoxon_direction_fidelity_yield_raw | 0.309235225 | 0.295875144 | 0.312449243 | +0.003214018 | yes |
| de_wilcoxon_direction_reach_raw | 0.156390888 | 0.153491703 | 0.155999025 | -0.000391863 | no |
| de_wilcoxon_sig_jaccard | 0.019730854 | 0.019882469 | 0.019646387 | -0.000084467 | no |

Raw expression error changes **-0.178%**, while raw fidelity yield changes **+1.039%**. The zero scaled expression score is a clipped score, not zero prediction error.

## What this experiment establishes

The V7 and V4 response-effect arrays match exactly. V7 selected the V4 transfer fallback because none of the expanded-source, pooling, neural or ridge alternatives passed the fixed calibration guards. Count generation uses the same zero-aware binomial/Poisson recipe with a new seed and different guide/depth-stratified control draws. This publication is a resampling comparison for the existing estimator; it does not measure the hidden performance of the rejected alternatives.

The overall difference is small. Two realizations cannot estimate a reliable seed-variation distribution or establish statistical significance. Aggregate scores do not identify the weak context or target. Keep V4 as the published reference and diagnose source-to-context response calibration and unsupported DE calls before adopting a more complex model. Require improvements across paired offline seeds; another unchanged-model upload is not a biological modeling advance.

## Provenance and reproduction

Submitted `2026-10-07T17:33:44.652819Z`; upload took 13 minutes 18 seconds. Rank 1050 at the recorded API check; ranks can move. Partition `val`, panel `vcc2026-val-1`, anchors `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`.

Generation-time readiness, model, package, independent-validation and pipeline records are preserved. Publication is stored separately in `leaderboard_result.json`, `published_api_status.json`, `leaderboard_comparison.csv`, `raw_metric_comparison.csv` and `publication_manifest.json`. The package SHA-256 identifies the independently checked local artifact; the status API does not expose the uploaded package hash.

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc status k0pkYuJRcw3C0U8K0NgL
.\.competition_env\Scripts\python.exe record_v7_result.py
```

These inspect the existing entry. [V7_SUBMISSION.md](V7_SUBMISSION.md) preserves the original model selection, validation and historical upload command.
