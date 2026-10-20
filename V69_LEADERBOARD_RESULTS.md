# V69 published result: count mixture rejected

Entry `dj3mrGx3Md59leld1ZI3` is **published**, according to the user-supplied CLI receipt. Overall **0.095262845** versus V66 **0.150882188**: **-0.055619343 (-36.86%)**. V66 remains the highest published project score. Rank 772 is a snapshot; ranks from different dates are not a controlled comparison.

| Scaled score (higher is better) | V66 | V69 | Change |
|---|---:|---:|---:|
| Overall | 0.150882 | 0.095263 | -0.055619 |
| PDS | 0.651701 | 0.472266 | -0.179435 |
| MSE | 0.096575 | 0.009662 | -0.086913 |
| NMAE | 0.085271 | 0.083253 | -0.002018 |
| FID | -0.022609 | unreported | unavailable |
| REACH | 0.079626 | unreported | unavailable |
| JAC | 0.014730 | 0.008060 | -0.006670 |

FID and REACH are absent from the supplied receipt, so their individual values remain null. MSE and NMAE above are scaled scores, not raw errors. The four reported component declines contribute **-0.045839102** to overall under the established six-metric equal-weight rule. PDS and expression accuracy account for most of the drop.

## What the experiment establishes

V69 changed 100 of 400 cell positions with a probability-coupled count bridge while retaining 300 exact V66 rows. Frozen V66 response means were generator inputs; this did not guarantee preservation of realized expression profiles. Official preparation and independent full-payload validity checks passed. A valid submission can still predict the biology poorly.

Public confirmation had already rejected the candidate: mean paired study profile-error ratio **1.107517** (about **10.75% worse**), with failed fidelity, reach and Jaccard guards and unavailable aggregate NMAE. Public PDS increased from **0.755556** to **0.822222**, but hosted PDS declined. Beating the matched shuffled mixture was insufficient evidence of a useful change. These reused public studies do not reliably rank candidates for hidden contexts.

This rejects this particular quarter-population bridge as an improvement over V66. It does not isolate one cause or disprove all distribution models. Mean calibration, sampling, source-assay mismatch and metric transfer are confounded. Another blend or a larger variance change is not supported by this result.

## Next decision

Retain V66 as the reference. Before another release, require improvement over V66 on held-out response profiles and discrimination, with adequate DE measurements and representative control nulls. Prioritize an assay-aware response-mean/count-likelihood teacher on an independent study or reserved donor; compare against V66 and shuffled evidence before generating cells. Freeze the split and thresholds before opening that reserve. This is a proposed experiment, not a completed model or a promised ranking gain.

## Reproduction and provenance

Run `python record_submission_result_69.py`, then `python record_submission_result_69.py --check` from the root. Only the standard library is required. JSON and CSV comparisons are in `analysis/v69_published/`; input/output hashes are in `publication_audit.json`. Exact reported decimal strings are retained.

The original scientific release remains sealed: 775 files, local package SHA256 `faf909db096d55d21e2fc1d0f194ea18dc94e3a8c42acd1c57673be9366e06b4`. The user's command links the uploaded artifact; the receipt supplies no server checksum. Historical pre-upload files still say not uploaded and remain unchanged; this publication record supersedes that status. The recorder performs no API query, training, prediction generation, upload or GitHub push.
