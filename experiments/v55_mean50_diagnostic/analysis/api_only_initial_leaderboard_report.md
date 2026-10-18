# V55 published result

Reviewed18 October2026 (Europe/Oslo). Official authenticated API confirms entry `CTdmP6ishgalaMMz4XRe`, model `Fixed V24 V36 half-blend diagnostic v55`, **published**, rank **575 at this check**. Historical ranks change with other submissions; compare scores on the same panel and anchors.

Overall **0.138767074**, change versus V24 **-0.004230781** (-2.96%), versus V36 **-0.003078408**. **V24 remains the best published project submission.**

| Scaled metric, higher is better | V24 | V36 | V55 | V55 minus V24 | Overall contribution vs V24 |
|---|---:|---:|---:|---:|---:|
| overall | 0.142997855 | 0.141845482 | 0.138767074 | -0.004230781 | n/a |
| pds | 0.644246367 | 0.643344247 | 0.643211404 | -0.001034963 | -0.000172494 |
| mse | 0.047568407 | 0.050173576 | 0.048972059 | +0.001403652 | +0.000233942 |
| nmae | 0.082270573 | 0.080934254 | 0.080649120 | -0.001621453 | -0.000270242 |
| fid | -0.015791837 | -0.016431778 | -0.016718192 | -0.000926354 | -0.000154392 |
| reach | 0.084839840 | 0.078160221 | 0.061690648 | -0.023149192 | -0.003858199 |
| jac | 0.014853782 | 0.014892370 | 0.014797405 | -0.000056377 | -0.000009396 |

REACH contributes -0.003858199 to overall, **91.2% of the net decline**. Expression accuracy improves, while the other five scaled components decline versus V24. This is arithmetic decomposition of the mean of six scores, not a causal attribution. MSE/NMAE columns above are scaled accuracies and increase when favorable; their raw errors decrease when favorable.

## What the diagnostic establishes

The fixed half-blend failed to improve on either parent. Moving halfway between response effects did not interpolate the hosted score between the parent scores. Count normalization, rounding and DE ranking/significance are nonlinear; this aggregate outcome does not identify which operation or genes caused the decline.

Public raw fidelity improved by +0.001319770, while private scaled fidelity changed by -0.000926354. The apparent4/4 local-to-official fidelity sign agreement from V54 now has a counterexample. It was a small, heterogeneous retrospective sample, not a reliable model-selection rule. Expression improved in both public and private comparisons, but this alone did not select a better overall model.

Public raw REACH changed by -0.001434685; private scaled REACH changed by -0.023149192. The direction of the loss was visible in the new public tests and in the original V37 rejection. Poor correlation of local REACH gains with private gains was not evidence that a REACH loss could safely be disregarded. Raw public/private scaled magnitudes cannot be equated.

Retain V24 as reference. This diagnostic supplies no evidence for another similar residual blend. Before choosing another upload, test the mean-response versus distribution/rank bottleneck on public data with fixed matched comparisons and independent reserve discipline. No new training, package, reupload or reserve opening is part of this result check.

## Provenance and preserved evidence

Partition `val`, panel `vcc2026-val-1`, anchors `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4` match V24 and V36. API submission timestamp: `2026-10-17T19:53:42.177382Z` (UTC). Local package SHA-256 `b877efdd31565ed7a58bcf8f26384a8837e546221fab5733ab0602e5db713c1b`; 4,073,697,280 bytes. The API does not expose an uploaded-file checksum; it confirms the owned model/entry, not remote byte identity.

The original [V55 generation report](V55_MEAN50_RESULTS.md), protocol, ready-file receipt, final manifest, arrays and package remain unchanged, including their generation-time unknown-score/uploaded-false fields. This publication record is separate. The [official API receipt](experiments/v55_mean50_diagnostic/analysis/published_api_status.json) and [canonical result](experiments/v55_mean50_diagnostic/analysis/leaderboard_result.json) record the subsequent outcome.

~~~bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc status CTdmP6ishgalaMMz4XRe
"C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe" record_v55_result.py --check
~~~
