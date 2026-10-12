# V16 published leaderboard result

V16 is published as `wG3GGbaTsl2dXFn7cQy1`, overall **0.141613857**, a gain of **+0.003730723** (2.71%) over V13 on the same panel and anchors. It is the best recorded project validation score. Rank was **562** in the user report and **561** at the API check; historical ranks are not directly comparable.

## Scaled comparison

| metric | v13 | v16 | change_from_v13 | overall_contribution_change |
| --- | --- | --- | --- | --- |
| overall | 0.137883 | 0.141614 | 0.003731 |  |
| pds | 0.639136 | 0.641553 | 0.002417 | 0.000403 |
| mse | 0.032878 | 0.046085 | 0.013207 | 0.002201 |
| nmae | 0.081822 | 0.080345 | -0.001477 | -0.000246 |
| fid | -0.015618 | -0.017296 | -0.001677 | -0.000280 |
| reach | 0.075574 | 0.084381 | 0.008807 | 0.001468 |
| jac | 0.013508 | 0.014615 | 0.001107 | 0.000185 |

All six displayed components are scaled scores: higher is better, including MSE and NMAE. Four improve; fidelity and fold-change accuracy decline. Expression accuracy contributes +0.002201 and reach +0.001468 to overall. Contributions are score changes divided by six, reconcile with the overall gain, and are arithmetic rather than causal attribution.

## Raw comparison

| metric | v13 | v16 | change_from_v13 | relative_change | higher_is_better | improved |
| --- | --- | --- | --- | --- | --- | --- |
| pds_cosine | 0.789487 | 0.790561 | 0.001074 | 0.001360 | True | True |
| expr_mse_unbiased_capped_norm | 0.976550 | 0.961835 | -0.014715 | -0.015068 | False | True |
| de_wilcoxon_lfc_nmae | 0.951212 | 0.952086 | 0.000874 | 0.000919 | False | False |
| de_wilcoxon_direction_fidelity_yield_raw | 0.508106 | 0.507558 | -0.000548 | -0.001079 | True | False |
| de_wilcoxon_direction_reach_raw | 0.145631 | 0.153714 | 0.008083 | 0.055507 | True | True |
| de_wilcoxon_sig_jaccard | 0.035250 | 0.035660 | 0.000410 | 0.011644 | True | True |

Raw expression error decreases 1.51%; raw reach increases 5.55%. Fidelity remains below its zero-score anchor. Improved count generation has not resolved response direction or biological dispersion.

## What we learned

Production retains the exact V13 response arrays and changes the population generator to residual10: 0.25 expected + 0.65 Poisson(expected) + 0.10 empirical response counts, followed by the original dual-moment solver. Gene-calibrated shrinkage and balanced molecule transport failed the frozen development guards and were not deployed. This result supports retaining residual10 as the new reference, with only modest evidence from one leaderboard experiment.

Local reserved comparisons improved raw NMAE by 0.00304047 on the eligible H1 panels, but the hosted NMAE score declines. Missing HCT116/HEK293T NMAE panels prevented complete six-metric verification. Fidelity also declines on the reserved and hosted comparisons. Local tests remain exploratory on partly reused biology; the later score does not retroactively turn them into independent or complete validation.

Residual10 retains only 20.5% of real total variance in the public control diagnostic and still produces substantial false DE there. Preserving control variance with balanced transport fixed that diagnostic but worsened actual perturbation accuracy. Next tests should hold residual10 fixed while separating response-direction errors from distribution errors with real treatment-replica controls, eligible per-family metrics, and matched seed/mean ablations. Another broad ensemble is not supported by these results: the V15 ensemble and V16 shrinkage alternatives were rejected.

## Provenance and reproduction

Submitted `2026-10-12T07:32:46.609932Z`. Upload took 16 minutes 26 seconds. Partition `val`, panel `vcc2026-val-1`, anchors `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`. Package: 4,073,615,360 bytes; SHA-256 `aec15ac5e5ef38e835bb4206f62c431772a28690fe5823490f045c5f0319cc08`. Official prep and independent full-payload/provenance/replay validation passed before upload. The API does not expose the uploaded file hash; the user's upload record connects the audited local artifact to this entry. This recorder checks file size and the prior sealed validation records, without rehashing the 3.79 GiB artifact.

The original report manifest, its 18 inputs, all 24 model dependencies, generation-time uploaded=false snapshots, and numerical code remain unchanged. Both original report outputs are archived byte-for-byte. `publication_manifest.json` maps their original hashes to those archives and separately records the current publication documents. After rerunning the generation-time report writer, rerun this recorder to restore the publication overlay. README is a project index and is not included in the immutable report hashes.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc status wG3GGbaTsl2dXFn7cQy1
.\.competition_env\Scripts\python.exe record_v16_result.py --api-file experiments/v16_control_preserving/analysis/published_api_status.json
```
