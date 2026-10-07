# V6 strategy: preserve response fidelity and calibrate cell distributions

This is an evidence-based proposal, not a generated submission or a predicted leaderboard score. V4 is the starting point because it is the best published candidate. V5 provides useful ideas and better expression/PDS performance, but its fidelity loss outweighs its other gains. The neural experiment supplies control-cell reconstruction models; it does not supply measured perturbation responses.

## Evidence that changes the next experiment

V5 improves five scaled components, yet loses 0.062963 overall against V4. Direction fidelity contributes -0.107541 to that difference; the other components contribute +0.044578. Raw expression error falls 22.25%, while raw fidelity yield falls 61.66%. Full records are in [V5_LEADERBOARD_RESULTS.md](V5_LEADERBOARD_RESULTS.md).

`analyze_v6_strategy.py` now compares all 900 frozen response groups on matching official axes. It excludes the own target from changes, and uses provided control mean CPM >5 for the direction comparison. The 0.05 absolute log2FC cutoff is descriptive, not a significance test.

| Context | Groups with changed mean effects | Changed downstream coordinates | Opposite signs among expressed effects >=0.05 in both | Median expressed-effect cosine |
|---|---:|---:|---:|---:|
| A | 22 | 22 | 0.00% | 1.000 |
| B | 300 | 3,084,007 | 0.77% | 0.905 |
| C | 300 | 3,083,712 | 0.83% | 0.906 |

A changes only 22 promoter-neighbor coordinates in 22 target groups. B/C change more broadly, but fewer than 1% of coordinates with sizeable effects in both versions reverse sign. Both versions' direction predictions can still be wrong relative to reality, and zero/new/weak effects are outside that percentage. These comparisons do not identify the cause of fidelity loss. They do make it especially important to isolate the generator, including distribution changes, normalized moments and significance calibration.

The official metric suite forms DE tables using Wilcoxon tests on per-cell CPM, a control-only CPM >5 filter and within-perturbation BH. Expression/PDS read pooled counts. Therefore matching group means and pooled profiles does not ensure appropriate DE calls. All six exclude the own target; PDS excludes all 300 panel genes. The six scaled components have equal weight. [Official scoring definitions](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics.md), checked 2026-10-07.

## Recommended candidate design

Start with V4's context-specific source policy: supported CD4 responses for A; calibrated direct H1 or transferred responses for B/C. Keep A/B/C baselines and sampling separate. Retain V4 as an unchanged comparator throughout. A/B/C cluster identity does not identify biological mixtures or establish causal source equivalence.

Use the V4 response vector as the default. Treat additional H1 measurements and V5 calibration as candidate updates, not universal replacements. Estimate reliability from independent donor/batch agreement and target-held-out source validation. For source-supported updates, evaluate modest effect-space blends, initially weights 0, 0.1 and 0.25 toward the alternative; these are a search grid, not selected settings. Conflicting directions should receive stronger shrinkage unless held-out source evidence supports the alternative. Do not eliminate all weak effects: fidelity also penalizes insufficient response coverage. Validate call budget, sign yield and set overlap together.

Separate response estimation from count emission. Begin with independently sampled real controls and the V4 thinning/addition generator. Evaluate a moment-calibrated alternative that targets moments implied by the SAME response vector. Do not graft V5's pooled deltas onto V4 mean effects: that changes both the response and the generator. Measure the emitted CPM fold changes, pooled profiles, zero/detection fractions, variance, expression-depth covariance and empirical rank distributions against the requested profiles.

Use neural models cautiously as cell-state/conditional-mean estimators. The selected raw architectures are residual in A and deep in B/C. Normalized neural improvements over ridge are modest; direct raw gains partly reflect the different objective. A proposed neural template mixes observed counts with conditional expected counts, with weight zero as the default and small nonzero weights considered only after validation. Recalibrate baseline means and pooled profiles and restore empirical residual variability rather than emitting smooth predictions as if they were observed cells. Learn/calibrate residuals on held-out control guides, not production-model in-sample errors. Rank features retain zeros and average ties.

Evaluate that neural template only after isolating the simpler generator. Prefer real controls if it creates baseline DE, reduces dispersion incorrectly, or worsens held-out perturbation direction yield. Use separate repeated guide splits/seeds for the template comparison and report that the previous neural test set has already been inspected.

Reducing one of the neural network's 300 inputs and propagating the result is not a validated intervention. The network learned observational relationships, including cell state and technical variation. Do not use that counterfactual as the default downstream response. A later neural response model would need training on measured perturbation/control contrasts and intervention-held-out validation.

## First experiments: isolate the components

| Experiment | Response | Generator | Purpose |
|---|---|---|---|
| Reference | Frozen V4 | V4 | Reproduce the strongest known candidate |
| Generator ablation | Frozen V4 | Dual moments derived from V4 response | Test distribution/moment fitting independently |
| Response ablation | V5 response | V4 | Test response replacement independently |
| Interaction comparator | V5 response | V5 | Existing published V5 combination |
| Neural extension | Best supported response | Selected generator plus small, calibrated neural template weight | Test added cell-state modeling after simpler choices |

Use identical original control draws, target/guide/depth balance and random seeds across ablations where possible. Evaluate more than one seed; repeated draws are uncertainty checks rather than extra independent biological datasets. Audit own-target knockdown and resulting library-normalization shifts, especially for abundant targets. The own target is excluded from scoring, but changing its counts can still change other genes' normalized expression.

## Validation that decides what is adopted

1. Split controls by guide group; assess no-response templates against disjoint observed controls. Measure mean/pooled drift, detection, variance/covariance and false DE beyond the unchanged-generator comparator.
2. Use actual public H1 perturbation/control cells with disjoint source control IDs/batches, and disjoint CD4 donors for effect validation. Hold all conditions of each target together in outer folds. Fit preprocessing, priors, blend strengths and emission settings inside training folds; evaluate once on outer folds.
3. Generate cells before judging DE. Use the official CPM/Wilcoxon/BH definitions, or the pinned scorer where compatible. Report all six raw metrics, direction-call budgets and actual emitted effects. Replace the old normal-SE approximation where raw source cells permit an actual rank test.
4. Use 400 distinct source cells where available. For smaller source groups, report reduced-size evidence separately; do not duplicate cells and treat them as independent observations to manufacture a 400-cell reference. CD4 summary releases do not support a complete cell-level generation benchmark by themselves.
5. Accept an update only if held-out source direction yield/overlap remain competitive and expression error improves consistently. Assess uncertainty across targets, source partitions and generation seeds, and retain the unchanged comparator when evidence is weak. Source-domain scores and self-consistency checks cannot guarantee a hidden A/B/C score.
6. Inspect weakest target/context groups before packaging. Verify all 900 groups, the exact 18,533-gene axis, 400 finite nonnegative integer-count cells per group, valid library sizes, reproducibility and official packaging. Freeze all settings and record hashes before leaderboard upload.

## Priority and fallback

The first V6 decision should be **V4 effects plus a separately validated generator**. This is a more interpretable test than replacing effects, templates and neural assumptions together. If generation improvements are unsupported, retain V4 generation and test only individually validated response updates. If neither improves source validation, preserve V4 and report the negative experiment rather than assuming a new package is better.

A wholesale V4/V5 cell blend can manufacture artificial heterogeneity; prefer effect-space updates followed by one calibrated generator. A larger control-trained network is a lower priority than a faithful perturbation benchmark.

## Reproduction and provenance

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe analyze_v6_strategy.py
```

The root diagnostic reads frozen V4/V5 response arrays and original context gene statistics, writes 900 target rows and three context summaries under `analysis/v6_strategy/`, and records exact model/code hashes and limitations. It does not alter previous models or generate/upload a package. These model comparisons complement, and do not replace, [V4_MODEL_STRATEGY.md](V4_MODEL_STRATEGY.md), [V5_MODEL_STRATEGY.md](V5_MODEL_STRATEGY.md) and [GENE300_NEURAL_RESULTS.md](GENE300_NEURAL_RESULTS.md).
