# V66 published result: small reach-driven gain

V66 published as `x8V3fUSgFIoAzHbNx4eX`, overall **0.150882188**, versus V63 **0.150422795**: **+0.000459393 (+0.31%)**. This is the highest published project score, but the gain is small and its statistical significance is unknown. Keep V63 as a near-tied comparison reference.

The user-reported rank is **518**, compared with V63's earlier 516. These are different-time leaderboard snapshots: other teams' results change, so a higher score can accompany a worse rank. Rank movement alone does not measure this model change.

## Component arithmetic

| Scaled metric | V63 | V66 | Change | Overall contribution |
| --- | ---: | ---: | ---: | ---: |
| pds | 0.651757503 | 0.651700951 | -0.000056552 | -0.000009425 |
| mse | 0.098470180 | 0.096574862 | -0.001895319 | -0.000315886 |
| nmae | 0.086290263 | 0.085270614 | -0.001019650 | -0.000169942 |
| fid | -0.018998265 | -0.022608595 | -0.003610330 | -0.000601722 |
| reach | 0.069881044 | 0.079625635 | +0.009744591 | +0.001624099 |
| jac | 0.015136044 | 0.014729659 | -0.000406385 | -0.000067731 |

All six CLI values are scaled, higher-is-better accuracies, including MSE and NMAE. Their arithmetic mean reproduces each overall. Reach is the only improvement. It contributes **+0.001624099** to overall; the other five metrics together contribute **-0.001164706**. Their sum is the small net gain. Fidelity is more negative, and expression, discrimination, fold-change accuracy and significance overlap all decline.

Against V24, overall gains **+0.007884332 (+5.51%)**: PDS, expression and fold-change accuracy improve; fidelity, reach and overlap worsen. The prior 18-publication review remains a frozen historical snapshot; V66 is the nineteenth publication now recorded.

## What changed and what the tests actually predicted

V66 raises the original assay-response proposal from half strength to three-quarter strength before applying the same V24 target-margin guard. The representation, source evidence and residual count generator are inherited from V63. The guard can reduce the accepted step, so the final effects are not uniformly multiplied by 1.5. See the sealed [method and results](V66_RESULTS.md). This whole-version comparison does not isolate a causal benefit from strength alone.

| Metric | Small public count test | Hosted validation |
| --- | --- | --- |
| pds | improves | worsens |
| mse | improves | worsens |
| nmae | improves | worsens |
| fid | improves | worsens |
| reach | worsens | improves |
| jac | worsens | worsens |

Five of six directions disagree; overlap alone agrees. Public error quantities are oriented so that smaller profile error and NMAE error count as improvements. Their magnitudes cannot be compared directly with hosted scaled accuracies. The public count test had 30 unique target/context pairs and two technical seeds; its 240 comparison records include shared baselines rather than 240 independent biological cases.

The paired profile-error ratio was **0.996334814** (about 0.37% lower error), concentrated in Orion; H1 and CD4 were worse. Leave-one-target/context-out sensitivity crossed a ratio of one. The predeclared strict REACH guard failed, yet hosted REACH is the sole gain. These observations show that passing most small-cohort proxies did not forecast the hosted component changes.

The broader 655 target/context-pair desired-mean diagnostic already warned of slightly worse profile error (equal-study ratio 1.000903825) and discrimination (PDS change -0.001063193). Those directions agree with the hosted expression and PDS declines. This broader check measures intended bulk means before cell generation, and it reuses public biology; it is not an independent official-score validation set.

## Next experimental priority

Treat V66 as a narrow validation gain, not a broadly better response model. Before another strength grid or ensemble, improve the selection benchmark: compare small and broad cohorts on exactly the same target/context pairs, replay intended means and generated counts separately, and test whether study-held-out component rankings transfer to a genuinely untouched public dataset. Freeze selection rules before observing that holdout. Retain shuffled-response and unchanged-generator controls. This result motivates those tests; they have not been run as part of recording this receipt.

No per-context/per-target scores, raw errors, confidence intervals or server-side package checksum were provided. Aggregate scores cannot identify which context or gene produced the gain, quantify significance, or establish performance on the final test panel.

## Provenance and reproduction

Source: the user's complete published CLI receipt; no additional API request. Partition `val`, panel `vcc2026-val-1`, anchors `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`, identical to V63 and V24. Submission timestamp is unknown. Upload duration was 16 minutes 18 seconds. Local package: 4,073,666,560 bytes, SHA-256 `be899434d81eac45bb8312077fe789088fe49b9b772971d589dcf5414e8f22f1`. The command links the local path to the entry; the server receipt supplies no checksum to independently confirm uploaded bytes.

The original 653-file V65/V66 scientific seal and prior 522-file V63 seal remain unchanged. Their pre-upload status/recommendation fields are historical. Separate publication records and the newest README block describe this later result. The publication audit hashes its source receipts, sealed diagnostics, recorder and three output files; README remains mutable.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe record_v66_result.py
.\.competition_env\Scripts\python.exe record_v66_result.py --check
.\.competition_env\Scripts\python.exe seal_v65_v66_release.py --check
```

Recording and checking do not change predictions, upload files or push Git.
