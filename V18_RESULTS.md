# V18 results: PIE did not earn a release

Completed 2026-10-13. **Both pretrained-expert blends failed the frozen development gates. No V18 submission was built or uploaded. The eight reserved target/family pairs remain unscored.** V16 remains the released reference.

The experiment completed **140 full-count population comparisons**: five arms on 14 target/family pairs and two seeds, with 400 distinct predicted/template cells and 400 genuine treatment cells per case. H1, HCT116 and HEK293T were held out from V16 response fitting in turn. Inference and scoring are separate; treatment values were never inputs to PIE.

## Family-balanced raw results

| Metric | V16 | Confidence PIE blend | 25% PIE blend | Better direction |
|---|---:|---:|---:|---|
| PDS | 0.777778 | 0.777778 | 0.652778 | Higher |
| Fidelity | 0.606476 | 0.604852 | 0.542392 | Higher |
| Reach | 0.503430 | 0.528948 | 0.383250 | Higher |
| Jaccard | 0.022365 | 0.021717 | 0.018625 | Higher |
| NMAE on eligible panels | 0.849526 | 0.863096 | 0.903545 | Lower |
| Paired expression numerator change | Reference | +6.61% | +650.30% | Lower |

These are public **raw metrics**, not the private anchor-scaled leaderboard scores. HCT116 had fewer than ten true DE genes, so NMAE is unavailable for both seeds and is not filled with zero. The expression change averages paired relative numerator changes within family then across families. The plain mean expression ratio is dominated by HCT116's tiny denominator and is not used as the headline result.

The confidence blend improved only reach. It missed the NMAE and expression guards and the requirement for at least three improved metrics. It beat neither amplitude-matched V16 nor target-shuffled PIE on two of fidelity/reach/Jaccard. The unconditional 25% blend worsened all six available criteria. Selection retained `v16` under the rule frozen **before scores**.

## What this taught us

A pretrained model can load correctly and still be an inferior response expert for this task. This particular checkpoint has no training-target overlap with our panel and no known identity text for A/B/C; V16 already transfers measured responses for known targets. On several H1/HCT116 cases PIE's raw covered-gene LFC errors were much larger than V16's. Its confidence head reduced the damage but did not add enough target-specific information beyond the negative controls. This is evidence about the tested checkpoint, input adaptation and blend, not a general claim that PIE cannot work.

Clipping and own-gene handling matter: the smoke example contained an approximately -26 log2 output, while the allowed blend clips to +/-3 and preserves V16 own-target handling. Published `p_de` scores are not calibrated competition FDR values. We did not tune a new confidence threshold using these outcomes or spend the reserved verification to rescue a failed candidate.

## Verification and reproduction

- Strict checkpoint load, supported FP32 GPU inference, finite outputs and bitwise repeat passed.
- Author missing-evidence/data-loader parity passed.
- All 79 archived author source Git blobs matched the pinned revision.
- Complete original NCBI token ranges for 299 available targets were reconstructed and independently checked; no token truncation was used.
- Every inference control-mean/axis/mask array was independently recomputed from original frozen NTC rows and matched exactly.
- Case identities, disjoint fit/reference controls, source-family exclusions and cache/metric hashes passed.
- All **28 V16 comparator count digests matched the previous independently tested V17 V16 comparator exactly**.

Evidence: `experiments/v18_pie_expert/analysis/development_summary.csv`, `development_panels.csv`, `development_scores.csv`, `development_expert_diagnostics.csv`, `selection.json`, `development_independent_audit.json`, `controls_input_independent_audit.json`, `asset_independent_audit.json`, and the frozen `plan.json`. Setup, asset retrieval and commands are in [V18_STRATEGY.md](V18_STRATEGY.md).

The next separately frozen experiment is V19 covariance-conditioned response transfer. It retains the failed V18 record and uses a different model rather than retuning PIE after failure.
