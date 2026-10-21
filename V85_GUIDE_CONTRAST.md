# V85: guide differences as a new source of response evidence

This experiment asks whether two different guides targeting the same gene can separate an on-target dose response from a response shared by both guides. It differs from V31, which averaged guide-specific treatment effects and retained between-guide uncertainty. V85 subtracts the guide responses, estimates their own-RNA dose difference, and predicts a response at 90% target suppression.

## Result

**The corrected estimator abstains on every target. No new prediction, count population or submission was produced.** This is a source-power rejection, not an observed leaderboard loss or proof that biological effects are absent.

The two existing D2 source extractions provide 19 eligible resting and 18 eligible stimulated targets. Eligibility requires exactly two guide IDs, at least two shared lanes, at least four cells per guide/lane, at least 20 NTCs per lane, and at least 15 cells per guide in each deterministic cell half. The initial metadata feasibility count was 19 per condition; applying the exact guide/lane and per-guide half-size requirements retains 18 stimulated targets.

The first implementation averaged within-lane log ratios. It reported 6 resting and 10 stimulated strong targets, but also **6 strong targets in each permuted-guide control**. This source-only null result exposed sparse log-ratio bias. No recipient prediction or response score had been computed. The original code, plan, tests, complete source fits and nulls remain preserved.

A documented amendment pools **linear within-lane CPM differences first**, then normalizes by the common matched NTC mean. This removes the logarithm and its offset from the dose/response contrasts. The corrected real and permuted sources both have **zero qualifying targets**. For example, resting SRSF5 has half-specific dose F statistics about 9.78 and 50.45; the fixed requirement is at least 10 in both halves. It was not lowered after seeing this result.

## Fixed estimator and assumptions

Split cells within guide/lane into equal, disjoint halves; omit any odd leftover. Use identical harmonic cell-count lane weights for dose and downstream response. Estimate each downstream contrast using the opposite half's own-RNA contrast, then reverse and average. Require the dose signs to agree, both absolute contrasts to exceed 15 percentage points of control expression, both F statistics to reach 10, and control own expression above 5 CPM.

The predicted fractional response is the two-sample dose ratio times -0.9, converted to log2 fold change. Approximate ratio uncertainty includes shared control normalization, worst-covariance bounds, split disagreement, the log derivative and a fixed variance floor. A fixed signal prior weights the correction. Own/panel outputs remain unchanged; missing measurements abstain.

This is an instrument-inspired predictive model. Guide-specific off-target effects, post-transduction selection, dose nonlinearity and expression measurement error can violate causal assumptions. The F screen is heuristic, not a calibrated weak-instrument test. The first-stage nulls were used to correct the source statistic, so the amended experiment must not be described as untouched preregistration.

[Inspre's primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12579002/) motivates using guide interventions as instruments. We did not implement its network estimator. The [RICE author repository](https://github.com/HowardGech/RICE) was inspected as related methodology; no RICE model was executed.

## Evidence and reproduction

- `v85_guide_contrast.py`: source fitting, weak-dose abstention, matched pooled/permuted controls and inference.
- `test_v85_contrast.py`: known response recovery despite a shared treatment component, weak/opposing-dose abstention, guide-swap invariance and disjoint halves.
- `audit_v85_experiment.py sources`: independently reconstructs selected raw-count statistics and checks every saved ratio estimate and source cell split.
- `report_v85_abstention.py`: records both source stages and verifies that all arms reproduce V66 bitwise on an H1 query and the D3 query. Since both source confidence arrays are zero, this abstention holds for every query, not only those two.

The final plan and source-only amendment are in `experiments/v85_guide_contrast/analysis/plan.json`; the original is in `engineering_history/log_stage_plan.json`. Original complete fits remain in `analysis/teachers`, corrected fits in `analysis/linear_teachers`. Full repeated cell-index traces and arrays are regenerable and ignored by Git; `source_summary.json` retains all target-level first stages, failures and hashes.

On a fresh checkout with the existing scientific dependencies and data, the final implementation can fit the corrected source bank directly: `test_v85_contrast.py`, `run_v85_experiment.py fit`, `audit_v85_experiment.py sources`, `report_v85_abstention.py`. Do not overwrite existing frozen fits. To replay the initial failed version, restore the three `log_stage_*py.txt` source snapshots and original plan in a separate checkout; their recorded hashes identify the exact code used. `run_v85_experiment.py amend_source` is the historical source-only transition and intentionally refuses to run after recipient evaluation.
