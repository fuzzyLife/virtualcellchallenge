# V87 results: conditional low-rank response operator

**No V87 submission package was generated.** The sole primary failed four of five frozen public mean checks. V66/V82 remains the best published prediction at overall 0.1508821875498449.

The BOLERO-inspired implementation was fitted and tested; upstream BOLERO inference itself was not run. See [BOLERO_REVIEW_V87.md](BOLERO_REVIEW_V87.md) for paper/code findings and direct-use requirements.

## Fixed public comparison

655 previously reused public target/context pairs; equal study votes for H1, Orion and CD4. Profile error is the squared error of the log1p(50,000-scaled) pseudobulk profile; PDS below is a local raw metric change. These are not six-metric anchored competition scores.

| Fixed output | Profile error / V66 | Error reduction | PDS change |
|---|---:|---:|---:|
| conditional25 | 0.998490009 | 0.1510% | -0.010995784 |
| constant25 | 0.998382771 | 0.1617% | -0.010545643 |
| conditional25_shuffled | 0.996633090 | 0.3367% | -0.008951509 |

| Study | Conditional error ratio | Constant error ratio | Shuffled error ratio |
|---|---:|---:|---:|
| H1 | 0.999436670 | 0.999181374 | 0.998808302 |
| Orion | 0.998669388 | 0.998709032 | 0.994255280 |
| CD4 | 0.997363970 | 0.997257908 | 0.996835687 |

## Prespecified checks

- profile: FAIL.
- worst_study: PASS.
- pds: FAIL.
- beats_constant: FAIL.
- beats_shuffle: FAIL.

Primary error improved only 0.1510%, below the 0.5% requirement. PDS fell 0.010996, beyond the 0.005 loss limit. Both constant conditioning and a fully refitted shuffled-response model had lower expression error. No other arm was promoted after these outcomes. The downstream donor/count gate was not opened; no new count populations, second seed, h5ad or vcc were produced.

## What this tells us

Source-target validation favored conditional over constant adapters for H1 and Orion model fits but not the CD4 fit. That did not translate into better recipient-study predictions. Smaller profile error with worse discrimination, especially when shuffled response identities reduce profile error more, is consistent with generic response regularization rather than improved target-specific transfer. It is not proof that every conditional model fails or that all biological signal is absent.

Across the five public contexts, the primary reduced between-target predicted response energy by 5.1-6.7%, while the shuffled model reduced it by 7.8-17.1%. These are predicted pseudobulk-response differences, not cell-to-cell biological variance. The reduction is consistent with contraction toward common responses; its desirability must be judged against actual prediction metrics. The full diagnostic is retained in response_diversity_diagnostic.json.

This test adds conditional weight generation across a frozen neural response backbone. It does not add the sequence regulatory information, compatible atlas encoder, or quantitative activity axis that powers the BOLERO results. We still have three source study groups per outer fold (three or four source cell-context families), with context and assay confounded. Complete-case output coverage was 5,922-6,507 of 18,233 nonpanel genes. The model predicts group means; it does not learn the within-population cell distribution.

All six frozen backbones selected epoch 100, the fixed training-budget boundary. No convergence claim is supported. Conditional adapter epochs were selected only on held source target identities, using a fixed 100-epoch budget. A later training-budget or feature revision would require a separate protocol; this negative pilot does not rule out the architecture.

## Verification and preservation

Numerical tests checked zero-B exact equality to the base network, explicit per-example collapsed-weight matrix equivalence, both factor gradients, absent base gradients and monotone-rank invariance. Six fitted models passed basis/target-partition/epoch-selection checks. Poisoning all quantitative response data from an excluded study left the entire allowed source bank exactly unchanged. All 15 public output arrays passed independent profile-margin checks. An independent metric audit replays per-target profile errors, explicit tied cosine-distance ranks, equal-study aggregation and the frozen decision. See the tracked JSON receipts for exact checks and hashes.

No earlier scientific plan, code or prediction was edited. Large fitted weights/arrays remain local and ignored; their hashes, source identities, full training histories, all public outputs metadata and all scores are versioned. V85/V86 were committed separately as 5f59360.

## Reproduction

Use the root scripts and scientific runtime commands in BOLERO_REVIEW_V87.md. After mean scoring, run `audit_v87_predictions.py predictions`, `audit_v87_predictions.py metrics`, and `report_v87_results.py`. Replays require the local scientific arrays/data whose hashes are recorded. Source inspection/download commands require public network access.
