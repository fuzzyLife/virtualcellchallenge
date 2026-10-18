# State-response research: V59, V60 and V61

Three separately frozen versions tested the requested state-abundance/within-state decomposition. V24 remains the best published project submission (overall 0.142997855). These are repeatedly examined public datasets; source cells and generation seeds are not independent biological studies. D3 remains unopened.

## What changed and why

- **V59:** four shared statistical RNA states on 512 stable, nonpanel genes, experimental abundance ratios and conditional responses. It reuses exact V52 disjoint raw views. Before scoring, CPM closure was corrected; before transfer, absent source genes were masked for pairwise-valid basal Spearman routing. Original code and both correction receipts are preserved.
- **V60:** the NTC-only state audit found CD4 coordinates strongly associated with library depth/detection. A separately frozen variant regressed nonpanel depth/detection covariates using NTCs only. This can remove biological size/complexity too; it is not proof of technical confounding. The representation change did not pass matched-control response tests.
- **V61:** pool matched state sufficient statistics across lanes before nonlinear estimation/shrinkage. The old tiny-lane oracle diagnostic averaged unstable small-population errors and is superseded by population-level reconstruction; originals remain available. V61 retains the exact V60 state bases, source exclusions, bounds, primary/shuffle definitions and count tests. Empty individual lane-state controls borrow their lane's overall control mean before pooled support checks; this imputation is explicit.

All three deployable predictors exclude the recipient's complete study from response fitting and representation selection. The independent auditor makes recipient treatment fields unavailable and requires bitwise-identical predictions. Matching target evidence in other studies is legitimate cross-context retrieval, not unseen-target generalization.

## Whole-study response results

| version | study | primary_over_V24_MSE | beats_shuffled_and_pooled | response_gate |
|---|---|---|---|---|
| 59 | CD4 | 0.990043 | False | False |
| 59 | H1 | 0.981287 | True | False |
| 59 | Orion | 0.995603 | False | False |
| 60 | CD4 | 0.989740 | False | False |
| 60 | H1 | 0.979072 | False | False |
| 60 | Orion | 0.995435 | False | False |
| 61 | CD4 | 1.007561 | False | False |
| 61 | H1 | 1.065896 | False | False |
| 61 | Orion | 1.002560 | False | False |

V59 completed 48 count comparisons on H1/HCT116 before stopping after its failed response gate; its remaining planned count cases were not run. V60 stopped after the failed response gate and has no count result claimed. V61 completed all 264 planned count comparisons to characterize the fixed release candidate, including the separate state-cell diagnostic. No partial pilot is represented as a completed 264-comparison experiment.

## NTC state-quality check

| source | PC | library_spearman_v59 | library_spearman_v60 | detected_genes_spearman_v59 | detected_genes_spearman_v60 |
|---|---|---|---|---|---|
| H1 | 1 | 0.341682 | -0.001977 | 0.262521 | -0.011909 |
| H1 | 2 | 0.469416 | -0.022100 | 0.437830 | -0.027424 |
| HCT116 | 1 | 0.284780 | -0.015622 | 0.206929 | -0.023532 |
| HCT116 | 2 | 0.635059 | 0.011760 | 0.688064 | 0.010143 |
| HEK293T | 1 | 0.094821 | 0.073699 | -0.025049 | 0.051738 |
| HEK293T | 2 | 0.606443 | 0.002331 | 0.657163 | 0.004308 |
| CD4_Rest_D2 | 1 | 0.871384 | 0.001136 | 0.888379 | 0.001886 |
| CD4_Rest_D2 | 2 | 0.906546 | 0.023387 | 0.929418 | 0.026335 |
| CD4_Stim48_D2 | 1 | 0.896057 | -0.009558 | 0.915605 | -0.005802 |
| CD4_Stim48_D2 | 2 | 0.918372 | 0.016305 | 0.942557 | 0.021122 |

The state figures use lane-centered control coordinates; visual overlap cannot establish biological equivalence. The population-oracle correction is an openly post-hoc diagnostic and never enters fitted release features or candidate selection. See `evaluate_state_population_diagnostic.py --version 60` for the corrected comparison on the original lane-wise estimator.

[Detailed V61 results](V61_RESULTS.md) · [Fixed release and upload instructions](V61_SUBMISSION.md) · [V59 protocol](V59_STATE_COMPOSITION_STRATEGY.md) · [V60 protocol](V60_DEPTH_ADJUSTED_STRATEGY.md) · [V61 protocol](V61_POOLED_STATE_STRATEGY.md).
