# V86 results: multivariate Bayesian response transfer

**Original qualification passed: False.** The fixed primary is `interval25`; no other arm is promoted after seeing outcomes. No upload package was generated. V66/V82 remains the best published project result at **0.1508821875498449**.

[Protocol and commands](V86_STRATEGY.md) | [Guide-contrast experiment](V85_GUIDE_CONTRAST.md)

## Different mechanisms actually tested

V85 used differences between guides as a source of on-target response evidence. Its corrected negative controls behaved properly, but every source target failed at least one fixed dose-evidence requirement, so it supplied exactly zero corrections. The original sparse-log failure, source-only correction and full first-stage summaries are preserved.

V86 learns a mixture of response-sharing covariances and computes a conditional posterior for an unseen context. It uses measurement masks and source uncertainties, excludes the entire recipient study, retains query-specific variation and tests a shared-only covariance comparator. This is an independent mash-inspired model, not a reproduction of mashr or a published competition solution.

## Synthetic and model checks

Synthetic held-context error / copying one source: **0.699546**. This validates the implementation under its simulated assumptions; it is not real-data accuracy.

Intervals use the conditional mixture mean and variance with a normal approximation; they are not exact mixture quantiles or established frequentist coverage. Mixture weights are plug-in estimates.

Six priors were fitted: genuine and independently target-shifted source evidence for H1, Orion and CD4 outer exclusions. Each uses 20,000 source target/output coordinates. These correlated coordinates are working-likelihood units, not independent biological observations.

| Outer exclusion | Shifted labels | EM iterations | Final maximum weight change |
| --- | --- | ---: | ---: |
| H1 | False | 100 | 0.00362166 |
| H1 | True | 100 | 0.00362674 |
| Orion | False | 100 | 0.00310088 |
| Orion | True | 100 | 0.00319360 |
| CD4 | False | 100 | 0.00115201 |
| CD4 | True | 100 | 0.00119297 |

All six fits reached the fixed 100-update limit before the 1e-7 weight-change stopping tolerance. This is a bounded pilot, not a fully converged fit. Complete source optimization is a remaining model-development step; no extra iterations or component settings were selected from recipient scores. Model metadata retain full weights and histories.

Independent checks replay Gaussian likelihoods against SciPy, conditional posteriors by direct solves, covariance positive semidefiniteness and mixture moments. Poisoning all excluded-study responses, masks, uncertainties, doses and profiles leaves every allowed source-bank field unchanged. All final prediction files also pass independent 80% identity-margin replay.

## All 655 reused public pairs

| Arm | Equal-study profile error / V66 | PDS change |
| --- | ---: | ---: |
| interval25 | 0.996265839 | -0.000536 |
| posterior25 | 0.989600007 | -0.011181 |
| shared25 | 0.995987093 | -0.000753 |
| interval25_shuffled | 0.993357747 | -0.009482 |

| Primary development check | Passed |
| --- | --- |
| profile | False |
| worst_study | True |
| pds | True |
| beats_shared | False |
| beats_shuffle | False |

H1, Orion and CD4 receive equal study weight; within Orion/CD4, contexts are averaged. All original source exclusions remain enforced, including every CD4 state and donor for CD4 queries. The biology and public outcomes have been reused, so this is development evidence. Complete context-level results are retained in `analysis/mean_results.json`.

## D3 final-count diagnostic

The original metadata-fixed 14 targets use 400 distinct actual cells, 400 distinct NTC templates, original template library totals and a separate 4,994-cell NTC reference. D3 was examined in V84 and is not blind. Current V86 fitting/conditioning excludes all CD4 responses; recipient treatment values never enter its predictor.

| Arm | Profile error | PDS | Fidelity | Reach | Jaccard | NMAE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| v66 | 0.006497 | 0.565934 | 0.518865 | 0.006449 | 0.152418 | 0.983722 |
| interval25 | 0.006478 | 0.565934 | 0.521370 | 0.008382 | 0.153311 | 0.981878 |
| shared25 | 0.006479 | 0.565934 | 0.516549 | 0.006999 | 0.152318 | 0.983218 |
| interval25_shuffled | 0.006467 | 0.560440 | 0.518889 | 0.007353 | 0.152471 | 0.983920 |

Primary profile error change: **-0.2907%**; PDS change **+0.000000**.

| Primary count check | Passed |
| --- | --- |
| profile_MSE | True |
| pds | True |
| fidelity | True |
| reach | True |
| jaccard | True |
| nmae | True |
| three_gains | True |
| beats_shared25 | True |
| beats_interval25_shuffled | False |

The fixed primary improves 5 of the six raw point estimates: profile_MSE, fidelity, reach, jaccard, nmae. Discrimination is unchanged. These gains are small: the shared-only model also improves expression, and the target-shuffled model has still lower expression error. The primary therefore fails the shuffled-expression requirement despite its favorable DE point estimates.


Raw public statistics above are not hosted anchor-scaled scores and are not averaged into an invented competition overall. These targets and this donor were previously used. A single technical seed is not independent biological replication.

Second seed executed: **False**. It is conditional on passing both original development and first count requirements. No failed thresholds were relaxed.

## Baseline reuse and numerical integrity

The test emitted **42 new populations** and reused **14 previously audited V66 populations**, for 56 compared records. Exact equality checks cover baseline response, CPM, bulk proportions, axes, cohort and seed; unchanged emitter/metric code is checked. Current truth/reference profiles are replayed before copying baseline evidence. The reused baseline is not claimed as newly generated data.

Emission checks 336 new raw tied-rank coordinates; 112 baseline coordinates were checked previously. Independent arithmetic replay verifies BH, directional scores, profile errors and 14-target discrimination for every one of the 56 records. Integer nonnegative counts and exact original row library totals are asserted.

## Decision and reproducibility

The fixed primary must pass development, count noninferiority and gain requirements, beat both controls, and pass a separately confirmed second seed before a release recommendation. Novelty and synthetic success do not satisfy those requirements. No new upload file or better hosted score is claimed.

Run `finish_v86_experiment.py` in a fresh analysis checkout with the committed plan and required prior data. Individual stages and numerical audits are in `V86_STRATEGY.md`. `report_v86_results.py` regenerates this report from the completed, hashed receipts. Existing scientific releases and their arrays remain unchanged.
