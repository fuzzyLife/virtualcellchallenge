# V30 full-target structured response results

**Decision: no V30 release.** V24 remains the best published project entry at 0.142997855. The full-target and construct-aware means failed the frozen screen; no generated counts, private scoring, fresh reserve or package followed.

## Across-context mean screen

|Variant|LFC MSE ratio to held-family V24|LFC MAE ratio|
|---|---:|---:|
|string_full|0.993255715|0.995970903|
|string_panel|0.993880035|0.996177439|
|random_full|0.992941008|0.995596610|
|construct_full|0.993255582|0.995972234|

Full STRING reduces squared error 0.67%, but randomized embeddings reduce it 0.71%. Both primary candidates miss the required 1% MSE and 0.5% MAE improvements and fail to beat the randomized control. They beat panel-only training by a very small amount. Construct-aware pooling changes the result negligibly at this precision. No statistical-equivalence claim is made.

|Family|Full STRING MSE ratio|Randomized MSE ratio|
|---|---:|---:|
|H1|0.991368433|0.991447564|
|Orion|0.998843517|0.998697087|
|CD4|0.989555195|0.988678371|

There are 31 active source/variant fits and one explicit RPE1 panel-only abstention; each active source fit has four rank/penalty solutions. The frozen nested selection records 96 inner family/configuration comparisons and 12 outer variant comparisons. These are summary-response tests on reused public biology, not independent biological replications or private anchor-scaled metrics. Whole-Orion and whole-CD4 exclusions are enforced.

## Separate unseen-target stress test

|Source|STRING MSE ratio to fitting-target mean|Randomized MSE ratio|Held targets|
|---|---:|---:|---:|
|K562_full|0.983436501|1.000705852|1648|
|RPE1_full|0.920977318|1.004673784|384|

STRING learns useful unseen-target signal within a measured source: squared-error reductions are 1.66% in K562 and 7.90% in RPE1 relative to the fitting-target mean. Randomized embeddings do not improve that baseline. This is a different task from predicting a new cell context. It cannot rescue the failed cross-context release candidate.

The target test was separately frozen after H1 context results were visible, before any target-stress outcomes. A target-identity hash creates 60/20/20 fitting/inner/test partitions; choose rank/penalty on inner only, refit 80%, then score 20% once. Fitting means, target sets and saved predictions are independently replayed. The secondary versus-zero diagnostic reproduces the original whole-vector float32 dot product; primary versus-constant errors use float64.

## Verification and scope

Independent original-row reconstructions match ten K562/RPE1 target summaries exactly. Five representative source/variant audits verify the masked-loss gradient, explicit abstention and saved predictions; full saved prediction replays have zero error. The worst checked relative masked-gradient residual is below 7e-9. Nested choices and actual excluded-family model requests pass. The separate target stress audit checks all four models, disjoint target roles, fitted bias and prediction arithmetic. Two solver fixtures pass, including equivalence to an independently constructed masked least-squares system.

This result applies to the implemented source-specific linear response bases, STRING features and rank-weighted transfer. It does not rule out nonlinear context-conditioned mean learning, conditional transport or foundation-model fine-tuning. Guide-resolved work is separately recorded in [V31_RESULTS.md](V31_RESULTS.md).

Reproduction: [V30_STRATEGY.md](V30_STRATEGY.md), root `collect_v30_embeddings.py`, `run_v30_experiment.py`, `run_v30_target_holdout.py`, `audit_v30_experiment.py`, `audit_v30_target_holdout.py`. The target stress outcomes are locked; use its auditor to verify an existing run rather than reopening selection. Large arrays remain outside Git.
