# Project retrospective: what worked, what failed, and what remains missing

Reviewed against local code and frozen evidence on 2026-10-11. The best published artifact remains V13, entry `AvAsYkuxFSjbxYDCd0wJ`, overall **0.137883134**. This review creates no model, prediction package or upload.

## Main assessment

We improved substantially, but the production model is still principally **a pooled public-intervention response plus a synthetic count generator**. The biggest useful changes were stronger measured response coverage (V4), then pooling/dose correction/post-pooling shrinkage with compact emission (V10), and consensus coverage/dispersion changes (V11). V13 corrects source knockdown handling and gains only **0.000879160** over V11. We are refining an effective but restrictive model class.

The most consequential new code finding is that the current model cannot learn context-dependent reversals in its pre-emission response vector. A/B/C controls are distinct, yet most expressed-gene responses are forced to be shared. The next major research step should test conditional intervention responses and identify the response-versus-distribution error budget. A larger control-reconstruction network does not resolve that by itself.

## Published results on the same panel

All ten recorded submissions use `vcc2026-val-1` and the same `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4` anchors. Compare scores rather than historical ranks, which change as others submit. The first record comes from user CLI output; subsequent records include read-only API confirmation. This retrospective made no API requests.

| Submission | Overall | Change from best before this submission | Main interpretation |
| --- | ---: | ---: | --- |
| First | -0.276290 | — | Conservative K562 transfer; weak downstream coverage |
| V2 | -0.223325 | +0.052964 | CD4 transfer and coverage helped |
| V4 | -0.050674 | +0.172652 | Lineage/source policy and stronger coverage helped |
| V5 | -0.113637 | -0.062963 | Fidelity loss outweighed five component gains |
| V6 | -0.056790 | -0.006116 | Local guards did not protect leaderboard fidelity |
| V7 | -0.051387 | -0.000713 | Same V4 effects, new sampling; no breakthrough |
| V9 | -0.070447 | -0.019773 | Same V4 effects, calibration lost fidelity |
| V10 | 0.116754 | +0.167428 | Pooled interventions plus compact populations |
| V11 | 0.137004 | +0.020250 | Consensus coverage plus dispersion tuning |
| V13 | 0.137883 | +0.000879 | Small KD-calibration improvement with tradeoffs |

V3 and V8 were generated diagnostic candidates without recorded leaderboard scores. V12 generated no submission. Version numbers do not count independent biological-model advances.

V10's discrimination/fidelity contributions explain **97.42% of its net gain over V4**. Response and emission changed together. V10/V11 have isolated local ablations, but the uploaded combined gains cannot be assigned to a single component without a matched leaderboard ablation. V11 versus V10 improves five scaled components while expression accuracy worsens. V13 versus V11 improves four while reach/overlap worsen. Its tiny gain has no repeated-submission uncertainty estimate and is not a statistically established breakthrough.

Exact scores, raw metrics where recorded, and additive component contributions are in [published_history.csv](analysis/project_retrospective/published_history.csv) and [component_attribution.csv](analysis/project_retrospective/component_attribution.csv).

## What we actually tested

| Work | What happened | What it establishes |
| --- | --- | --- |
| Full/300-gene rank analysis, mapping and markers | A/B/C separate strongly; broader markers support T-like A and epithelial-associated B/C; external units are displaced from all three clouds | Baseline context structure is real. Nearest-reference expression is not a validated causal-transfer weight or exact cell-line identity. |
| V1/V2/V4 coverage and lineage transfer | Large sequential gains, with multiple changes bundled | Measured interventions and correctly directed downstream coverage matter. More calls alone are not enough. |
| V3 basal deconvolution/ensembles | Generated, unsubmitted | Expression-fit improvement never demonstrated leaderboard improvement; coefficients are not cell-type fractions. |
| V5/V9 moment and population calibration | Five scaled components improved versus V4, but fidelity worsened enough to lower overall | Better mean/bulk calibration can still hurt significant-call yield. |
| 300-to-other-gene ridge and neural models | Neural raw-objective MSE skill about 26–32% against depth baseline; normalized gains over ridge modest | Control reconstruction demonstrates association. Reducing a measured input is not validated CRISPRi inference. The raw-network gain changes the loss as well as architecture. |
| Tahoe/PCA response heads and V7 supervised neural/ridge alternatives | Particular interventional heads/blends failed generated-count guards | These small heads were not useful. Full STATE/CellOT and a proteomics-trained encoder were not trained/tested; foundation-model approaches are not generally ruled out. |
| V7 fallback and V8 evaluator repair | V7 response arrays equal V4; old PDS could reward constant predictions | Resampling is not a new causal model. Early near-perfect local discrimination was invalid. |
| V10 response pooling and population transport | First positive overall score, large discrimination/fidelity gain | A strong interventional baseline emerged; population bias remained. Integer moment transport is not learned neural optimal transport. |
| V11 consensus rescue plus noise 0.65 | Better coverage/fold changes/overall; expression worsened | Useful bias/variance tradeoff, without establishing biological single-cell realism. |
| V12 distribution variants | Empirical/gamma/response-linked/thinning/noise methods improve null realism but fail perturbation fidelity/error guards | Uniformly restoring control variance is not the missing recipe. Effects, detectability and treated heterogeneity must be evaluated together. |
| V12 verified K562 Flex and platform harmonization | Six source variants failed; pathway screen did not beat shuffled annotations | Matching platform or adding these annotations alone did not solve transfer. Flex is another K562 assay, not another independent lineage. |
| V13 low-rank programs and robust/KD variants | Low-rank finalists failed; Student3 promising but rejected; unknown-KD correction selected | The released V13 is not a response-program network. Some inexpensive alternatives remain unconfirmed. |

The full scope/evidence inventory is [experiment_inventory.csv](analysis/project_retrospective/experiment_inventory.csv). Relevant detailed reports include [GENE300_NEURAL_RESULTS.md](GENE300_NEURAL_RESULTS.md), [EMBEDDING_RESPONSE_PILOT.md](EMBEDDING_RESPONSE_PILOT.md), [V7_MODEL_STRATEGY.md](V7_MODEL_STRATEGY.md), [V12_RESEARCH.md](V12_RESEARCH.md) and [V13_STRATEGY.md](V13_STRATEGY.md).

RPE1 was parsed and audited: zero overlap with these 300 intervention targets, despite shared output genes. It cannot simply provide direct responses for this panel. Reactome and RNA-fingerprinting code were collected; fingerprint fitting and responder modeling were not executed. K562 construct-level KD metadata was recovered after V13's plan froze, so its quantitative values are not used by the released model. More source data are not universally useless, but overlap, assay and independent context diversity must be checked before promising improvement.

## Where the current score is weak

These are scaled scores and their contributions to the six-metric mean, not raw accuracies.

| Component | V13 scaled score | Contribution to overall |
| --- | ---: | ---: |
| Discrimination | 0.639136 | 0.106523 |
| Expression accuracy | 0.032878 | 0.005480 |
| DE fold-change accuracy | 0.081822 | 0.013637 |
| Direction fidelity | -0.015618 | -0.002603 |
| Direction reach | 0.075574 | 0.012596 |
| Significance overlap | 0.013508 | 0.002251 |

Discrimination supplies about **77.3% of the net overall score**. Its raw value is 0.78949, while raw significant-set Jaccard is only **0.03525**, raw direction reach **0.14563**, and corrected expression error **0.97655**. We distinguish some target profiles better than we recover their exact downstream DE sets. Raw fidelity 0.50811 is a response-budget-normalized yield, not 50.8% sign precision; Jaccard is not recall. These aggregate results cannot tell us which A/B/C context or targets fail.

V13 reduces raw expression error about 1.08% versus V11 while raw reach falls about 2.59%. Improving expression alone is not enough to improve every biological metric. See [V13_LEADERBOARD_RESULTS.md](V13_LEADERBOARD_RESULTS.md).

## Gap 1: context interactions are almost absent from production

[fit_v13_submission.py](fit_v13_submission.py) builds one all-family response bank. [v11_methods.py](v11_methods.py) calls the expression gate in [v8_methods.py](v8_methods.py). The query-context adjustment is:

```text
effect[context, target, gene] = pooled_effect[target, gene]
                              * clip(control_CPM[context, gene] / 5, 0, 1)
                              * clip(control_CPM[context, target] / 5, 0, 1)
```

The own-target downstream entry is zero; direct knockdown is emitted separately. Nonnegative gates can attenuate or silence effects, but cannot reverse their sign. There is no learned target-by-context interaction in this production equation and no context-specific response amplitude above the gate.

The new saved-array audit finds **288/300 targets** and **7,991/18,533 genes** have fully open gates in all three contexts. All **2,301,120 downstream target/gene pairs** in that subset have exactly identical A/B/C effect values. Every context pair has **zero opposite-sign nonzero effects**, across approximately 3.8 million joint-nonzero coordinates. Removing the expression gates recovers the same underlying vector to maximum absolute discrepancy **3.14e-7**, consistent with float32 rounding.

This is a representational limitation, not proof that particular hidden genes reverse direction. Different baselines, compositional normalization and synthetic draws can still yield different final fold changes or Wilcoxon calls. V4 used different lineage policies, and V7 tested simple context-conditioned heads. Neither supplies an accepted, validated interaction learner in V13.

Evidence: [context_effect_geometry.json](analysis/project_retrospective/context_effect_geometry.json), [context_effect_pairs.csv](analysis/project_retrospective/context_effect_pairs.csv), [target_context_gates.csv](analysis/project_retrospective/target_context_gates.csv).

## Gap 2: the current emitter does not learn cell states or responder fractions

`v11_methods.emit` starts from context-wide mean and pooled proportions, uses sampled real library depths, forms `expected = depth * desired_mean_probability`, mixes Poisson noise at amplitude 0.65, then integer-fits both moments. [build_v13_submission.py](build_v13_submission.py) records genuine control IDs and guide balance, but does not use their individual expression patterns as gene-expression templates. It retains depth provenance, rather than the actual sampled biological states.

The moment fitter creates dependencies and depth/profile coupling, so the final genes are not literally independent. There is nevertheless no learned treated-state covariance, responder fraction, cell-cycle redistribution or target-dependent residual population model. All groups share the same noise setting. Earlier real-template and thinning/immigration methods were tested and often underperformed; simply copying controls again is insufficient.

Existing zero-response diagnostics measured approximately 0.23 times real-control total CPM variance and about 1,711 null DE calls against 2,000 controls for V11. Empirical/gamma alternatives made controls more realistic but often destroyed perturbation fidelity. This is an unresolved tradeoff, not evidence that null realism guarantees a better leaderboard score. Source-summary LFC accuracy is therefore insufficient to predict count-level metrics. See [POST_V11_RESEARCH.md](POST_V11_RESEARCH.md) and [V12_RESEARCH.md](V12_RESEARCH.md).

## Gap 3: source confidence and coverage remain heuristic

Five families differ in assay, timing, cell counts and measurement masks. Fixed family weights and scalar disagreement penalties do not provide a calibrated target/context directional posterior. The post-V11 ungated summary audit found fold-change MSE worse than zero effects in all five held families, despite some transferred directions. This normal-q diagnostic is not the official generated-count score.

The same audit reports a median **4,677 non-target outputs without usable response support per target**, including abundance/QC exclusions as well as unmeasured genes. Removing masks indiscriminately could increase false calls. We need measured-versus-missing distinctions, uncertainty calibrated across held contexts, and explicit abstention or justified extrapolation. V13 fixes unknown KD being mistaken for zero but does not solve context-dependent sharing. See `analysis/post_v11_research/transfer_family_summary.csv` and `uncertainty_bins.csv`.

## Gap 4: validation arithmetic improved, but biology and power remain limited

Historical issues were found and repaired: target-specific PDS origins; subset-only fidelity protection; a V10 validation-depth adapter; missing official NMAE eligibility; and exact Scanpy tie/rank arithmetic. The current CPU scorer has pinned official-function and independent count fixtures. Old results are historical evidence, not equally valid tests of current metric semantics. See [TRANSFER_VALIDATION_AUDIT.md](TRANSFER_VALIDATION_AUDIT.md), [EVALUATION_CRITERIA_AUDIT.md](EVALUATION_CRITERIA_AUDIT.md), and the corrected-depth section of [V10_STRATEGY.md](V10_STRATEGY.md).

That does not remove the transfer gap. Only three source families have generated-count comparisons; CD4/K562 evidence is mainly summaries. Technical strata and seeds are not independent lineages. Many real treatment groups contain roughly 150–200 cells while predictions contain 400, so DE power differs from the challenge. Private A/B/C treated labels and matched anchors are unavailable: public overall proxies are not official validation scores.

V13 development used **58 target/seed cases**, including reused H1, and confirmation **30** cases from **HCT116/HEK293T only**. Some confirmation panels had no NMAE-eligible contrasts, making the complete primary guard inconclusive. Original confirmation outcomes have now been viewed and cannot be called unopened. The expanded V12 plan contains **8 supplemental confirmation target/seed identities** absent from V13; this review compares identities without evaluating those counts or establishing adequate power. Larger control samples or duplicated treatment cells cannot create independent treatment information.

Repeated development on the same panels also risks adapting to public-source peculiarities. The next benchmark needs larger genuine treatment groups where available, held-family/assay splits, training-only transformations, and a newly frozen final check whose identities/power are established before candidate comparisons. Summary responses previously inspected are not pristine blind biology. See [validation_scope.json](analysis/project_retrospective/validation_scope.json).

## Gap 5: our selection rule can reject the competition's preferred tradeoff

V13 required strict expression-error improvement, then chose the passing candidate with the greatest expression gain. The competition equally averages six context-scaled metrics. These objectives differ.

The Student-t source-reliability variant improved **all five other primary development metrics**, including reach **+0.022280**, while expression ratio worsened **0.002999, approximately 0.213%**. Its first-seed public-source anchored composite was **0.030403 versus 0.013840** for V11. It was rejected correctly under the frozen rule, but was not demonstrated inferior under every reasonable objective. It has no independent count confirmation or leaderboard result. Reused panels and small public-anchor denominators make it an unconfirmed lead.

Do not retroactively relax guards and rename development outcomes confirmation. Prospectively compare a bounded Pareto set using all six metrics, retain raw-metric and worst-family safeguards, and predeclare a defensible expression noninferiority margin from replicate/sampling uncertainty. This review does not establish an arbitrary acceptable margin or promise that Student3 will beat V13. Evidence: [student3_policy_tradeoff.json](analysis/project_retrospective/student3_policy_tradeoff.json) and [selected_local_ablation_history.csv](analysis/project_retrospective/selected_local_ablation_history.csv).

## Next experiments should answer different questions

These are proposed tests, not executed results or a new submission. Keep exact V13 and V11 comparators. Isolate changes before combining winners.

| Order | Experiment | Question and useful decision |
| --- | --- | --- |
| 1 | Response-versus-population oracle diagnostic | In held public contexts, compare transferred versus independently estimated true treated mean/bulk responses under the same compact emitter; separately test a treatment-derived residual-distribution oracle. If correct means still lose DE metrics, generation limits us; if means rescue most loss, prioritize transfer. Treatment information in an oracle is diagnostic only and isolated from deployable fitting/selection. Estimate and score on disjoint treatment cells where possible; a same-truth oracle is explicitly optimistic. |
| 2 | Conserved response plus context interaction | Fit `conserved[target,gene] + interaction(target, baseline_programs, gene)` using actual intervention contrasts. Start with a regularized interaction model against pooled V13, zero, generic-response, target-shuffled and matched-capacity context-free controls. Learn program mappings inside held-family folds. Promote only if target-specific residual transfer and six count metrics improve in fresh contexts. |
| 3 | State-dependent population response | Learn treatment-linked state reweighting and within-state response/responder variation from source cells; map query controls into training-defined programs. Compare composition changes alone, within-state changes alone, and their combination with fixed response means. Test NTC nulls and genuine held treatments. This differs from global gamma/Poisson noise or blindly copying controls; absent query treatments cannot identify their fingerprints directly. |
| 4 | Strongest inexpensive unconfirmed head | Freeze Student3 against V13/V11; test recovered construct KD separately, then combine only if supported. Use prospectively defined objective/tolerances and new powered count confirmation. Do not select from the spent V13 confirmation or treat five development gains as proof. |
| 5 | Data acquisition with identifiable benefit | Inventory independent CRISPRi contexts with target overlap, usable axes and genuine cell power. Prioritize context diversity and comparable source controls. Broader intervention dictionaries may train interaction programs without exact target overlap, but need a separate target-held extrapolation test. More K562 platforms alone do not supply unseen-context evidence. |

For tests 1–3, distinguish generic stress/viability responses from target-specific information using generic-centroid, shuffled-target and NTC pseudo-intervention controls. Preserve official common-control scoring and report residual diagnostics separately. Covariance/state realism is explanatory evidence until intervention tests demonstrate benefit; the leaderboard does not directly score every covariance.

The practical priority is **a powered error-budget experiment, then conditional response learning**. Preserve Student3 for a fair, inexpensive prospective comparison. Another broad architecture or dispersion sweep on the same panels is lower priority.

## Reproduction and audit limits

Run from Windows Command Prompt in the project root:

```cmd
.\.competition_env\Scripts\python.exe review_project_experiments.py
```

Root [review_project_experiments.py](review_project_experiments.py) checks all ten score means and panel identities, exports exact component changes, compares stored V11/V13 ablations, verifies V13 model-array and frozen-dependency hashes, and audits every saved A/B/C effect coordinate. It compares validation-plan identities without running old evaluators or reopening count confirmation. It records source/output hashes and actual runtime locations, including existing user-site scientific dependencies.

The report's narrative is human-reviewed; rerunning the script regenerates its supporting CSV/JSON evidence and records the report's hash. Restore the ignored V13 model array from the frozen generation pipeline if absent; a hash mismatch fails the audit. Inputs are checked again after the audit to detect concurrent changes. The manifest records exact coverage in [audit_manifest.json](analysis/project_retrospective/audit_manifest.json).

This verifies mathematical and metadata findings, not hidden biological accuracy or fresh package integrity. It makes no claim about other teams' current undisclosed methods. Prior literature and participant reports remain evidence in the linked project reports; no new internet survey was performed for this retrospective. All frozen models, plans, packages and publication records are preserved.
