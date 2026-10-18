# V61 published leaderboard result

V61 published as `JNIsuxo0KEVQveGuZoLX`, overall **0.116763474**, compared with V24 **0.142997855**: **-0.026234381 (-18.35%)**. V24 remains the recommended published reference. Rank **674** is the user's snapshot; historical ranks change as other teams submit.

## Same-panel comparison

| Metric | V24 | V61 | Change | Contribution to overall change |
| --- | ---: | ---: | ---: | ---: |
| pds | 0.644246367 | 0.556224522 | -0.088021845 | -0.014670307 |
| mse | 0.047568407 | 0.000000000 | -0.047568407 | -0.007928068 |
| nmae | 0.082270573 | 0.082669197 | +0.000398624 | +0.000066437 |
| fid | -0.015791837 | -0.018906992 | -0.003115155 | -0.000519192 |
| reach | 0.084839840 | 0.069947269 | -0.014892571 | -0.002482095 |
| jac | 0.014853782 | 0.010646848 | -0.004206933 | -0.000701156 |

All six components above are scaled, higher-is-better scores, including MSE and NMAE. Their equal-weight average reproduces each reported overall; component changes divided by six reproduce the overall decline. The V61 expression score of zero means it did not earn a positive scaled expression improvement over the baseline. It does not mean zero raw prediction error. Raw V61 metrics were not supplied and are not inferred.

PDS and expression together account for **86.1%** of the net decline. NMAE is the sole gain, only +0.000398624; fidelity, reach and overlap also worsen. This arithmetic locates the lost score, but does not isolate a causal mechanism or identify context A, B or C as responsible.

## Interpretation and next decision

The state-response approach improved within-dataset reconstruction for CD4 and H1 but failed its held-study response and count-generation support gates before upload. Its equal-study response MSE ratio was 1.025339 versus V24 (lower is better), and count-profile expression error worsened in every held study. The hosted regression agrees with that warning. See [frozen public results](V61_RESULTS.md) and [V59-V61 review](V59_V61_STATE_REVIEW.md).

The evidence does not support stronger weighting of this state correction or an ensemble that assumes it adds independent useful signal. Retain V24 as the reference. Future candidates should pass the predeclared held-study response and integer-count gates before another submission is recommended. These reused public datasets cannot establish final-test generalization, and this one aggregate leaderboard result does not prove the exact cause of failure.

## Provenance and reproduction

Source: the user's complete published CLI receipt, without an additional API query. Partition `val`, panel `vcc2026-val-1`, anchors `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`, matching the saved V24 receipt. No submission date was supplied. The recorded local package hash is `cf4807ba59156619ae1bc2a315a87450bab5a1021c5c8cc569703b059118ad53` (4,068,689,920 bytes); the user's command links the entry to this path, but the receipt supplies no server-side artifact hash. Upload duration: 16 minutes 27 seconds.

The original scientific code, model, package, pre-upload reports and 297-file release manifest are preserved. Their historical uploaded=false fields describe generation time. The separate publication receipt, comparison CSV and publication audit record this subsequent outcome. README is the current status overview.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe record_v61_result.py
.\.competition_env\Scripts\python.exe record_v61_result.py --check
.\.competition_env\Scripts\python.exe seal_v61_release.py --check
```

Reproduction writes only publication outputs and the marked README status block. It does not fit models, alter predictions, submit a package or push Git.
