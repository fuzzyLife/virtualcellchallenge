# V14 experiment results (2026-10-11)

**Decision: reject the context-interaction candidate; retain V13.** This research experiment does not create or upload a prediction package. V13 remains the best recorded project submission, overall 0.137883134. Public raw scores below are not private anchor-scaled leaderboard scores.

## What was tested

V14 retains the exact held-family V13 pooled response and adds a bounded, target/gene-specific interaction with two control-expression rank components. The context feature basis and all response slopes exclude the entire evaluated family. CD4's conditions remain one family. No panel target genes enter the feature basis, missing measurements abstain, and own-target effects remain excluded. Penalties 2 and 8 were screened; penalty 2 advanced. Fixed shuffled context labels and simple 1.1x amplification serve as controls. See [V14_STRATEGY.md](V14_STRATEGY.md) and root `v14_context_methods.py`.

Actual-count development compares four models across nine targets in H1, HCT116 and HEK293T, two paired seeds, 400 distinct real treatments and 400 predictions per case. Full independent controls are primary; sampled400 references are secondary. That is **72 synthetic populations**, 144 scored population/reference combinations, and 48 profile-panel rows. Seeds are repeated draws of the same nine biological contrasts, not independent lineages. Models use public training families only; held-family treatment data enter evaluation, never fitting.

Separately, nine oracle cases split 400 distinct treatment cells into 200 calibration and 200 truth cells. Three arms produce **27 population cases: 18 synthetic populations and nine independent real-treatment replicas**. Combined with development, there are 99 population cases, including 90 generated populations. Oracle calibration/truth splits are disjoint within each case, but reuse development targets and can overlap development identities across experiments; they are not independent confirmation experiments.

## Response screen

The best interaction proxy barely exceeds V13. Both shuffled labels and amplification score higher, which does not support interpreting the proxy gain as learned context biology. K562's held-out feature coverage is only 41.5%, below the predeclared 70% threshold; it abstains and cannot demonstrate an interaction benefit.

| candidate | mse_ratio | fidelity_proxy | jaccard_proxy | sign_reference_DE | selection_objective |
| --- | --- | --- | --- | --- | --- |
| amplitude110 | 1.13994 | 0.301774 | 0.0644752 | 0.611325 | 0.141941 |
| context2 | 1.12674 | 0.288185 | 0.0653167 | 0.605825 | 0.136577 |
| context8 | 1.12459 | 0.285662 | 0.0658776 | 0.611225 | 0.136084 |
| generic | 1.00181 | 0 | 0.170702 | 0.527551 | 0.0422232 |
| shuffled2 | 1.13161 | 0.289985 | 0.0664733 | 0.602776 | 0.137067 |
| target_shuffle | 1.15862 | 0.241592 | 0.0264299 | 0.435542 | 0.099277 |
| v13 | 1.12504 | 0.285545 | 0.0662862 | 0.611325 | 0.136182 |
| zero | 1 | 0 | 0.170702 | 0 | 0.0426755 |

These are previously inspected public summary contrasts, not a blind test or official count metrics. Zero/generic controls can have misleadingly high set-overlap proxies for weak/empty call sets.

## Completed 400-cell count test

Higher is better for PDS, fidelity, reach and Jaccard; lower is better for NMAE and expression ratio. Scores first average within source family, then across families with available values. **NMAE is unavailable for all HCT116 primary panels**, so its displayed aggregate covers H1/HEK293T only. Missing primary metrics block promotion rather than being filled with zero.

| candidate | pds | fidelity | reach | jaccard | nmae | common_expression_ratio |
| --- | --- | --- | --- | --- | --- | --- |
| amplitude110 | 0.833333 | 0.608503 | 0.439115 | 0.0110986 | 0.88434 | 197.149 |
| context2 | 0.722222 | 0.605056 | 0.455305 | 0.0111558 | 0.853533 | 201.563 |
| shuffled2 | 0.833333 | 0.602542 | 0.452158 | 0.0107345 | 0.848803 | 196.564 |
| v13 | 0.833333 | 0.610206 | 0.395159 | 0.0113139 | 0.851627 | 194.604 |

| metric | context2_minus_v13 |
| --- | --- |
| pds | -0.111111 |
| fidelity | -0.00514996 |
| reach | 0.0601461 |
| jaccard | -0.000158163 |
| nmae | 0.0019056 |
| common_expression_ratio | 6.9589 |

Only reach improves. The candidate fails these frozen checks: **primary_complete, expression_noninferiority, pds, at_least_four_metrics_improve, utility_improves, beats_shuffled_context_utility**. Primary discrimination drops 0.111111; fidelity drops 0.005150; Jaccard drops 0.000158; NMAE increases 0.001906. The expression ratio worsens beyond the prospectively declared exploratory 1% tolerance. The unanchored utility also loses to V13 and shuffled labels. This utility is dominated by HCT116's tiny expression denominator and must not be read as an estimated competition score. Independent failures in discrimination, metric breadth and missing-metric completeness support retaining V13 without relying on that utility.

### Family results

| family | candidate | pds | fidelity | reach | jaccard | nmae | common_expression_ratio |
| --- | --- | --- | --- | --- | --- | --- | --- |
| H1 | context2 | 1 | 0.65488 | 0.126777 | 0.0232859 | 0.872629 | 1.44986 |
| H1 | v13 | 1 | 0.659585 | 0.0574493 | 0.0240342 | 0.875911 | 1.43091 |
| HCT116 | context2 | 0.5 | 0.616534 | 0.805556 | 0.000521514 | NA | 602.416 |
| HCT116 | v13 | 0.5 | 0.618711 | 0.694444 | 0.000380737 | NA | 581.666 |
| HEK293T | context2 | 0.666667 | 0.543753 | 0.433584 | 0.00965981 | 0.834437 | 0.823942 |
| HEK293T | v13 | 1 | 0.552321 | 0.433584 | 0.00952687 | 0.827343 | 0.716004 |

### Metric power and denominator diagnosis

The baseline's mean real DE-call counts, eligible fraction of the three target cases, and NMAE for each primary panel are:

| family | seed | reference_calls | nmae_eligible | nmae |
| --- | --- | --- | --- | --- |
| H1 | 2.0261e+09 | 156.667 | 1 | 0.877014 |
| H1 | 2.0261e+09 | 164.333 | 1 | 0.874808 |
| HCT116 | 2.0261e+09 | 2.66667 | 0 | NA |
| HCT116 | 2.0261e+09 | 2.66667 | 0 | NA |
| HEK293T | 2.0261e+09 | 29.3333 | 0.666667 | 0.843625 |
| HEK293T | 2.0261e+09 | 29.3333 | 0.666667 | 0.811061 |

Expression numerator and denominator are retained separately. HCT116 has an extremely small control-distance denominator; enormous ratios are unstable and must not be used as a headline biological error ranking. Family ratios are calculated on common-origin panels rather than averaged per-target ratios.

| family | candidate | common_expression_ratio | common_expression_numerator | common_expression_denominator |
| --- | --- | --- | --- | --- |
| H1 | context2 | 1.44986 | 0.00167749 | 0.00115693 |
| H1 | v13 | 1.43091 | 0.00165552 | 0.00115693 |
| HCT116 | context2 | 602.416 | 0.00104213 | 1.72992e-06 |
| HCT116 | v13 | 581.666 | 0.00100623 | 1.72992e-06 |
| HEK293T | context2 | 0.823942 | 0.00130406 | 0.0015827 |
| HEK293T | v13 | 0.716004 | 0.00113322 | 0.0015827 |

No reserved V14 algorithm-verification outcomes were opened after rejection. Earlier history already contains some of those biological contrasts, so they cannot be described as blind biological confirmation in any future run.

## Treatment-mean versus population diagnostic

The oracle uses **known treated-calibration depths in every arm**, including the V13 comparator. This deliberately isolates response/population limitations; it is not deployable baseline inference. Measured treatment mean-CPM and pooled proportions constrain compact emission; actual replicas use the original 200 calibration cells. All are evaluated against the disjoint 200 truth cells using the same independent controls and scorer.

| family | candidate | nmae | jaccard | predicted_calls | reference_calls |
| --- | --- | --- | --- | --- | --- |
| H1 | v13_compact | 0.890997 | 0.0151542 | 319.333 | 61 |
| H1 | oracle_moments_compact | 0.240116 | 0.0260288 | 1643.33 | 61 |
| H1 | actual_treatment_replica | 0.236306 | 0.296327 | 60.3333 | 61 |
| HCT116 | v13_compact | NA | 0 | 263.667 | 1 |
| HCT116 | oracle_moments_compact | NA | 0.000724341 | 1453 | 1 |
| HCT116 | actual_treatment_replica | NA | 0.666667 | 0.666667 | 1 |
| HEK293T | v13_compact | 0.569759 | 0.0440119 | 137.667 | 12.3333 |
| HEK293T | oracle_moments_compact | 0.166933 | 0.00760392 | 1230.67 | 12.3333 |
| HEK293T | actual_treatment_replica | 0.162097 | 0.669516 | 12.3333 | 12.3333 |

For H1, supplying measured treatment moments reduces NMAE from **0.8910 to 0.2401**, close to the independent real replica's **0.2363**. But compact oracle cells call roughly **1,643 DE genes** against **61** mean truth calls, with Jaccard **0.0260**. Independent real treatment replicas call **60.3**, with Jaccard **0.2963**. The means help response accuracy while the synthetic distribution still fails to reproduce the DE set.

HEK293T shows the same distribution warning: oracle NMAE improves from **0.5698 to 0.1669** on its eligible contrast, but its mean DE calls increase to **1,231** versus **12.3** in truth. Jaccard falls from **0.0440 to 0.00760**, compared with **0.6695** for real replicas. Two of its three truth cases do not meet NMAE's ten-DE-gene minimum. HCT116 averages only one truth DE gene at 200 cells and has no eligible NMAE case; its ratios and discrete overlap estimates have little power.

This supports **two limitations**, not a proof that distribution alone explains our leaderboard plateau. Transferred means remain imperfect. Compact populations also behave very differently from independent real replicas even when given measured moments. DE calls beyond the truth set are not all proven false positives: the 200-cell truth has limited power, moments themselves are estimated from 200 cells, and residual batch composition can matter. Fitting moments does not fix zero frequencies, heterogeneous conditional means or tails. Expression error includes jackknife terms, so measured moments do not guarantee a perfect expression ratio.

## What this changes about the next experiment

Do not further tune this two-component interaction against the same few families. Five assay-confounded families are weak supervision for per-target conditional effects, and the shuffled-label control provides no persuasive evidence that this version learns useful biology. Additional independent interventional contexts, rather than more observational A/B/C control reconstruction, would improve identifiability.

Prioritize a **distribution diagnostic before another submission**. Use cross-fitted real treatment marginals conditioned on library-depth/state strata as a diagnostic oracle; preserve the joint ranks of empirical control cells while changing only source-supported per-gene quantiles. Compare that with the same means/pooled proportions and compact emitter, including a no-perturbation arm and shuffled target assignments. Test both 200 and 400 genuine cells. Require calibrated control-null DE calls, acceptable expression numerators, and improved fidelity/overlap across adequately powered held-family treatments. Oracle-treated distributions remain diagnostic, not legal inputs for hidden competition treatments.

Only if that diagnostic closes the oracle gap should we fit a deployable cross-context transport of **treated-versus-control quantile changes**, with gene-specific dropout/variance changes, context support checks and shrinkage to V13. This is more specific than globally restoring variance: uniform gamma/full-variance approaches already failed V12. Keep means, distribution changes and source selection as separate frozen ablations; validate before packaging. These are proposed next tests, not completed improvements or a promised leaderboard gain.

## Checks, provenance and limitations

Four scientific tests passed: a known context sign reversal, invalid/missing measurements not voting zero, poor-coverage abstention, and held-family label poisoning not changing predictions. The independent audit passed **27 cached profile files**, all planned original cell identities and disjoint splits, and **2592 independently recomputed ridge coordinates**. It replays scalar weighted normal equations, verifies coefficient centers/own-target exclusion, confirms family-balanced aggregations/case counts, and verifies frozen V13 numerical dependencies remain unchanged. It does not independently rerun every DE test or establish hidden biological correctness.

The original plan's freshness inventory missed V7 `cell_split_*.json` history. The separate `metadata_correction.json`, produced before numerical selection/count outcomes were inspected, supersedes that annotation. No identities, seeds, numerical methods or candidates changed. Eight of nine development target names were previously count-tested. Of four reserved target names, only MYRF was globally unused; CAT, DCAF16 and PSMB9 were not. Original freezes are retained to preserve the audit trail.

Controls are disjoint and split within batches, but treatment truth pools source batches at family level. Assay/lineage and residual batch composition remain confounded. Public small panels are noisy; no official private anchor set is available. No new biological dataset was acquired for V14. The 1% expression tolerance is an exploratory guard, not statistical noninferiority.

### Reproduction and files

Run the experiment commands in [V14_STRATEGY.md](V14_STRATEGY.md). With completed hash-verified caches, assessment is preserved once written; it refuses overwriting. Then run:

```cmd
.\.competition_env\Scripts\python.exe audit_v14_experiment.py
.\.competition_env\Scripts\python.exe document_v14_experiment.py
```

Root scripts include preparation, interaction fitting, execution, assessment, history audit, independent audit, documentation and focused tests. `dispatch_v14_oracle.py` only schedules the unchanged oracle evaluator after each model cache finishes; `oracle_dispatch.json` pins this wrapper. Large model/profile caches are ignored by Git. Plans, source/model sidecars, scores, selection, audit and report hashes are retained under `experiments/v14_context_interactions/analysis/`. `report_manifest.json` pins this document's generator and input/output bytes.
