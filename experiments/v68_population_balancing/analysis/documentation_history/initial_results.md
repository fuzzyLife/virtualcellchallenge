# V68 results: preserve real cells, test population selection

The fixed real-cell population model **did not qualify for release**. Paired equal-study generated-profile error / V63 is **1.213016** (+21.30%). Its raw direction fidelity is **0.004661**, versus V63 **0.574376**. The complete fitted null guard passes: **True**. No V68 submission was generated. V66 remains the highest published project score, 0.150882188.

## What was tested

This is an alternative population generator, not another response blend. The optimizer selects weights on individual real query-control cells to approach V63's desired bulk expression, with capped inclusion probabilities and fixed depth-group quotas. Pivotal sampling selects distinct original cells. Every non-own-gene count remains intact; only the own gene receives an independent 0.1 binomial thinning before fitting/sampling. There is no synthetic Poisson population, mean refit or per-gene count repair.

The desired means are shared with V63; **realized means are not forced to match**. Controls may not contain the states required to express the predicted treatment response. This experiment tests the whole alternative decoder, not a pure variance change with equal realized means. See the [frozen method and limitations](V68_STRATEGY.md).

## Treatment outcomes

| Raw diagnostic | V63 | Balanced real | Shuffled weights | Uniform real |
| --- | ---: | ---: | ---: | ---: |
| pds | 0.833333 | 0.500000 | 0.666667 | 0.500000 |
| fidelity | 0.574376 | 0.004661 | 0.000000 | 0.000000 |
| reach | 0.194507 | 0.094430 | 0.120956 | 0.159877 |
| jaccard | 0.062476 | 0.004439 | 0.000000 | 0.000000 |
| nmae | 0.935484 | 0.946650 | 0.952672 | 0.942458 |
| profile_square_error | 71.977551 | 90.285886 | 88.026649 | 93.239416 |

PDS, fidelity, reach and overlap are raw higher-better diagnostics. NMAE and profile-square error are lower-better errors. The paired profile ratio averages three study ratios; it differs from the ratio of the displayed mean absolute errors. Undefined case metrics are omitted within their context, and a metric remains undefined if an entire study is unresolved. These values are not the leaderboard's scaled accuracies, and no overall competition score is inferred. Ten target/context pairs are tested across five contexts and three studies. Four arms produce 40 comparison records, including ten hash-verified reused V63 outcomes. Actual treatment populations/references are identical across arms. Predicted groups have 200 or400 cells, matching the saved comparator; resting CD4 uses200.

Per-study profile ratios: H1 1.028114, Orion 1.352834, CD4 1.258099. Fixed checks failed: profile, pds, fidelity, reach, jaccard, nmae, two_DE_gains, beats_shuffled_weights_profile. No parameter, control arm or mixture was promoted after outcomes. Two-target discrimination panels are coarse and unstable; their ranks cannot establish a hosted gain.

In H1, the two fixed targets HDAC8 and NFE2L1 have 157 and137 reference DE calls; balanced-real predictions make zero downstream calls in each. Restoring natural baseline variation did not reproduce their measured induced transcriptional responses. This is direct evidence for these tested cases, not proof that all occupancy models fail.

## Control-null calibration

| Context | V63 null calls | Uniform-real null | Complete fitted null |
| --- | ---: | ---: | ---: |
| H1 | 1198 | 0 | 0 |
| HCT116 | 5328 | 0 | 0 |
| HEK293T | 545 | 0 | 0 |
| CD4_Rest_D1 | 964 | 0 | 0 |
| CD4_Stim48_D2 | 2775 | 0 | 0 |

Uniform sampling tests the raw real-cell sampler. The additional fitted zero-response test exercises the complete balancing pipeline. The [separate addendum](V68_NULL_AUDIT.md) was frozen after some treatment outcomes but before any null outcome, because review found that uniform sampling alone did not test optimization. It adds a requirement without altering the original candidate or guards. These are fifteen null-comparison records across five contexts, including the original ten and five additional fitted-null records; they are not fifteen independent studies.

![Population and null diagnostics](experiments/v68_population_balancing/analysis/population_tradeoff.png)

## What we learned and what remains open

Technical calibration and biological response are separate problems. Selecting intact baseline cells can preserve authentic dispersion, but cannot invent perturbed states absent from the fit-control population. This test does not support replacing V63's induced-response generator with the tested occupancy-only model. It also does not justify a hybrid mixture or stronger weight fit without a new protocol. The relevant next distinction is between baseline variability and variation created by a perturbation; broadening natural variation alone is insufficient in these cases. This experiment does not isolate whether the loss comes from the control-state convex hull, fixed regularization, imperfect upstream response means or finite-sample power; an oracle diagnostic would be needed to separate them.

V67 also failed: its separate network mechanism reduced broad mean-profile error by only0.0018% and lost the target confirmation check. Together, these experiments provide two explicit negative results rather than an unsupported new upload. They do not establish which method will improve the hidden panel. D3 remains unopened; the public biology here has been used before.

## Verification and reproduction

Eight V68 numerical tests pass, including6,000 sampling simulations, objective finite-difference checks and recovery of a known population shift. All five control frames match frozen V63 query/control bulk profiles. Control/reference/treatment partitions are disjoint. All ten weight fits are feasible; **4/10** meet the fixed1e-8 dual-gap stopping rule. Maximum final dual gap is **1.86e-06**. Effective sample sizes range **672.6 to 4371.1**, always at least the number of selected cells. The remaining fits reached the preset iteration limit, so this trial is not evidence that an exactly optimized occupancy model must fail. The reported dual gaps bound suboptimality of the convex surrogate, not uncertainty in the discrete DE scores. No solver or regularizer was changed after outcomes.

The independent audit reconstructs all30 new treatment populations from original control rows, verifies non-own count preservation, selected-row uniqueness and depth quotas, independently reproduces full profiles and count-derived fold changes, checks BH adjustments, and replays four DE metrics including raw-p-value tie ordering for reach. It manually verifies16 Wilcoxon rank coordinates per population. Cached V63 baseline receipts remain bound to their original hashes; old raw-p reach ties are not independently replayed from caches that did not save p-values. Null call totals are reconciled; no independent full regeneration of the old synthetic-null generator is claimed.

Root commands: [V68_STRATEGY.md](V68_STRATEGY.md), [additional null commands](V68_NULL_AUDIT.md), `audit_v68_results.py`, `document_v68_results.py`, and `seal_v68_experiment.py --check`. The original protocol was committed as `78fb103`; the stricter pre-null addendum as `96f5665`. Scientific scripts and compact records are committed locally; large arrays are ignored and regenerable. Prior scientific release seals remain unchanged. No upload or Git push occurs automatically.
