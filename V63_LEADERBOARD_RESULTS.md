# V63 published leaderboard result

V63 published as `mhKv7sajwD0DKQ8oHPYP`, overall **0.150422795**, compared with V24 **0.142997855**: **+0.007424940 (+5.19%)**. It is now the best published project reference, superseding V24 on this validation panel. User-reported rank **516** is a historical snapshot. V61 scored 0.116763474; V63 gains +0.033659321 over it.

## Same-panel comparison

| Metric | V24 | V63 | Change | Contribution to overall change |
| --- | ---: | ---: | ---: | ---: |
| pds | 0.644246367 | 0.651757503 | +0.007511137 | +0.001251856 |
| mse | 0.047568407 | 0.098470180 | +0.050901773 | +0.008483629 |
| nmae | 0.082270573 | 0.086290263 | +0.004019691 | +0.000669948 |
| fid | -0.015791837 | -0.018998265 | -0.003206428 | -0.000534405 |
| reach | 0.084839840 | 0.069881044 | -0.014958797 | -0.002493133 |
| jac | 0.014853782 | 0.015136044 | +0.000282263 | +0.000047044 |

All six values are scaled, higher-is-better accuracies, including MSE and NMAE. Their equal-weight averages reproduce the reported overalls; each change divided by six gives its contribution to the overall gain. Four metrics improve and two worsen. Raw V63 errors were not supplied: the higher expression score does not establish that raw error halved.

Expression contributes **+0.008483629** to overall, exceeding the net gain because fidelity and reach losses offset part of the improvements. PDS, log-fold-change accuracy and significance overlap also improve. Fidelity remains negative, and reach falls -0.014958797; these remain clear gaps relative to V24.

## What this supports

V63 changes the source-response estimator: it removes source-common response components, translates effects through an assay-scale expression representation, and limits the correction using V24's predicted target margins. It retains V24's count generator. The published result supports retaining this complete candidate as our current validation reference. It does not isolate the benefit of any individual component; that requires component ablations with comparable validation.

The public paired profile error ratio was 0.968979699 (3.10% lower plug-in error), consistent with the direction of the hosted expression improvement. However, public count discrimination worsened while hosted PDS improved. Reach improved in development but worsened in additional public cases and in this submission. The strict public gates failed before upload and remain failed. This hosted win does not retroactively validate them: they are imperfect selection proxies, with limited unique targets and reused datasets. See the frozen [public results](V63_RESULTS.md) and [validation audit](V62_V64_VALIDATION_AUDIT.md).

## Next decision

Use V63 as the published comparison baseline, keeping V24 as the prior reference. Preserve the expression improvement while testing fidelity/reach corrections; require paired, study-held-out component ablations and shuffled controls before attributing a gain to a mechanism. Check predicted means and generated counts separately, since the public expression and discrimination checks disagree in their ability to forecast the hosted result. More complex ensembles are not justified by this result alone.

The receipt provides no per-context or per-target breakdown, so it cannot identify which context or gene caused the improvement. One aggregate validation outcome provides no estimate of uncertainty or guarantee for the final test panel.

## Provenance and reproduction

Source: the user's complete published CLI receipt, without an additional API query. Partition `val`, panel `vcc2026-val-1`, anchors `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`, matching both saved comparison receipts. Submission time was not supplied. Local package SHA-256: `a075397d8c23fa2e63e1e99220a04bf6fcd500eaa1e94f8a0fc8eb445271048b`; 4,073,656,320 bytes; upload duration 16 minutes 24 seconds. The user's command links this path to the entry, but the server receipt supplies no artifact hash.

The original scientific scripts, model, package, pre-upload reports and 522-file release manifest remain unchanged. Their uploaded=false and recommendation fields record pre-upload knowledge. Separate publication files record this subsequent result. README's newest status block is the current overview; older sections are historical.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe record_v63_result.py
.\.competition_env\Scripts\python.exe record_v63_result.py --check
.\.competition_env\Scripts\python.exe seal_v63_release.py --check
```

The recorder writes only publication files and the newest marked README block. It does not change predictions, upload files or push Git. Its check mode verifies recorded arithmetic and hashes; the release-seal check separately verifies the original scientific files.
