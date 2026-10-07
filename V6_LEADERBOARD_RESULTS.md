# V6 published leaderboard result

Entry `6MCjBWpPoCt5Sbw5CMu3` is **published**, overall **-0.056790**, rank 1052 at the user-reported check. V4 remains the best project submission at **-0.050674**. All five published project entries use the same validation panel and anchor set.

| Scaled metric | First | V2 | V4 | V5 | V6 | V6 minus V4 | Contribution to overall change |
|---|---:|---:|---:|---:|---:|---:|---:|
| overall | -0.276290 | -0.223325 | -0.050674 | -0.113637 | -0.056790 | -0.006116 |  |
| pds | 0.019644 | 0.103042 | 0.255302 | 0.447856 | 0.259374 | +0.004072 | +0.000679 |
| mse | 0.000000 | 0.000000 | 0.000000 | 0.035294 | 0.000000 | +0.000000 | +0.000000 |
| nmae | 0.016006 | 0.032352 | 0.067539 | 0.070189 | 0.076131 | +0.008592 | +0.001432 |
| fid | -1.665331 | -1.470413 | -0.683855 | -1.329100 | -0.730255 | -0.046400 | -0.007733 |
| reach | 0.033874 | 0.041380 | 0.086849 | 0.121061 | 0.083531 | -0.003318 | -0.000553 |
| jac | -0.061932 | -0.046314 | -0.029878 | -0.027120 | -0.029522 | +0.000356 | +0.000059 |

Overall changes by **-0.006116** versus V4 and **+0.056847** versus V5. Three scaled components improve (PDS, LFC accuracy and overlap), two worsen (fidelity and reach), and expression remains at its clipped zero score. Fidelity contributes **-0.007733**; the other five together contribute **+0.001617**. This is arithmetic attribution within the six-component mean, not causal attribution.

| Raw metric | V4 | V5 | V6 | V6 minus V4 | Improved versus V4 |
|---|---:|---:|---:|---:|---|
| pds_cosine | 0.616462 | 0.702367 | 0.618079 | +0.001616 | yes |
| expr_mse_unbiased_capped_norm | 1.247733 | 0.970079 | 1.253228 | +0.005495 | no |
| de_wilcoxon_lfc_nmae | 0.959164 | 0.957501 | 0.954593 | -0.004572 | yes |
| de_wilcoxon_direction_fidelity_yield_raw | 0.309235 | 0.118548 | 0.295875 | -0.013360 | no |
| de_wilcoxon_direction_reach_raw | 0.156391 | 0.186364 | 0.153492 | -0.002899 | no |
| de_wilcoxon_sig_jaccard | 0.019731 | 0.020722 | 0.019882 | +0.000152 | yes |

Raw expression error rises **0.44%**; direction-fidelity yield falls **4.32%**. A zero scaled expression score does not mean zero prediction error. Fidelity yield uses a response-call/reference budget; it is not simply the fraction of correct signs.

## What this changes

V6 retained the V4 stochastic generator and all A mean effects, blending 25% toward V5 only at existing same-sign B/C coordinates. Neural templates and full dual-moment generation were excluded. The positive public H1 direction/overlap results did not translate to a better aggregate anonymous-context leaderboard score. The source benchmark changes an H1 contrast gain, whereas production combines lineage estimators and applies an additional sign/zero mask; these are not identical response policies. Nine primary H1 targets and previously inspected source batches cannot establish anonymous-context transfer validity.

The control draw/emission seed also changed across all contexts. There is no published unchanged-V4 comparator using that seed, no context-level outcome breakdown, and no repeated-seed leaderboard uncertainty estimate. The aggregate difference cannot identify which context/target failed, establish statistical significance, or isolate the response blend from sampling variation. This result is not a test of neural perturbation prediction because the production neural weight was zero.

Keep V4 as the production reference. Future work should strengthen source-to-context validation and evaluate target/context-specific response reliability and seed stability before adopting another broad blend. V6 recovers much of V5's lost fidelity, but has not surpassed V4.

## Publication and provenance

Submitted at `2026-10-07T10:08:17.032649Z`. Upload took 13 minutes 29 seconds. Rank 1052 at the recorded API check; ranks change as other teams submit. Partition `val`, panel `vcc2026-val-1`, anchors `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`.

The generation-time `submission_ready.json`, frozen model, predictions and original numerical validation records are preserved. Publication is recorded separately in `leaderboard_result.json`, `published_api_status.json`, `leaderboard_comparison.csv` and `raw_metric_comparison.csv`. The package SHA-256 comes from the verified local readiness record; the status API does not independently expose or confirm the uploaded hash.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc status 6MCjBWpPoCt5Sbw5CMu3
.\.competition_env\Scripts\python.exe record_v6_result.py
```

These inspect the existing entry. [V6_SUBMISSION.md](V6_SUBMISSION.md) preserves generation, validation and the historical upload command.
