# V15 results: weighted response and population transport

No new package was released; V13 remains the upload reference. The completed comparisons did not justify replacing it under the frozen release policy.

V13 remains the best *published* project score, 0.137883134. Public raw metrics below do not predict an official leaderboard score or a rank gain. Numerical V13/V14 dependencies and prior packages were preserved.

The compact response blend changed fidelity by -0.001270, reach by -0.003969, Jaccard by -0.000032, and eligible NMAE by -0.002245. Its paired expression error worsened by 2.23%. These are family-balanced public changes; lower NMAE is better.

The tested ensemble did not supply enough independent error correction. The distribution variants increased apparent DE calls in control-only comparisons, and their response-count improvements did not survive the full guard policy. No reserved verification or production fit is opened after a failed development gate.

## Implemented experiment

Three response experts were tested: exact V13, Student-t pooling with neutral unknown-knockdown handling, and a training-source-reliability/similarity specialist. A six-point convex grid retains >=50% V13, with weights calibrated through inner-family holdouts inside each outer-held-family evaluation. This is weighted continuous-effect prediction, not a majority vote on counts. Missing measurements abstain; CD4 conditions remain one family. Per-gene direction agreement is diagnostic evidence, not independent model votes.

Source-treated/control marginal changes were extracted from verified H1/HCT116/HEK293T counts. Five relative quantiles in two depth bins are transferred onto ranked empirical query-control cells, then blended at 10% or 25% on measured supported genes and integer-fitted to the same response moments. This tests one gene-specific response-distribution construction. The linear blend also attenuates the original Poisson term by (1 - blend): its coefficient falls from 0.65 to 0.585 or 0.4875 on supported genes before fitting. It therefore changes existing sampling noise as well as adding empirical variation. The zero-response null uses zero quantile deltas and still fails, so these results do not isolate a failure of cross-source quantile transfer itself. A future noise-preserving construction would need a new frozen experiment, not a retrospective change to V15. Extreme quantiles, pooled source batches, assay differences and finite treatment counts remain limitations. See [V15_STRATEGY.md](V15_STRATEGY.md) for the exact equations, thresholds, identities and commands.

Population-case counts: **{'null': 18, 'development': 168}**. These include repeated seeds, multiple predictions for shared truths, and null/control diagnostics; they are not that many independent biological experiments. Development uses 14 targets with 400 genuine cells per case; reserved algorithm verification uses nine other current target names. Older public biology was already inspected, so this is not blind biological confirmation.

## Ensemble screen and complementarity

Frozen response for count testing: **student50**.

| candidate | mse_ratio | fidelity_proxy | jaccard_proxy | sign_reference_DE | selection_objective |
| --- | --- | --- | --- | --- | --- |
| ensemble | 1.08955 | 0.295967 | 0.0664901 | 0.607323 | 0.142426 |
| specialist | 1.06647 | 0.257437 | 0.0686852 | 0.61108 | 0.129385 |
| student | 1.10789 | 0.306913 | 0.064554 | 0.607541 | 0.143164 |
| student50 | 1.08955 | 0.295967 | 0.0664901 | 0.607323 | 0.142426 |
| target_shuffle | 1.12388 | 0.245201 | 0.0334579 | 0.447699 | 0.101917 |
| v13 | 1.0798 | 0.278998 | 0.0664159 | 0.608227 | 0.136248 |

Nested mixture weights:

| held_family | V13_weight | Student_weight | specialist_weight |
| --- | --- | --- | --- |
| CD4 | 0.5 | 0.5 | 0 |
| H1 | 0.5 | 0.5 | 0 |
| HCT116 | 0.5 | 0.5 | 0 |
| HEK293T | 0.5 | 0.5 | 0 |
| K562 | 0.5 | 0.5 | 0 |

Error correlations and direction corrections/breaks are averaged per source family/expert. The direction counts describe reference-DE genes per target/context contrast, not all-gene accuracy or official fidelity. High residual correlation limits the diversification benefit of an ensemble.

| family | expert | error_correlation | corrects_v13_sign | breaks_v13_sign |
| --- | --- | --- | --- | --- |
| CD4 | ensemble | 0.997298 | 0.344123 | 0.407155 |
| CD4 | specialist | 0.9994 | 2.05281 | 2.13799 |
| CD4 | student | 0.990136 | 0.643952 | 0.781942 |
| H1 | ensemble | 0.9978 | 5.63636 | 6.04545 |
| H1 | specialist | 0.990093 | 76.5455 | 65.9545 |
| H1 | student | 0.991575 | 10.8636 | 10 |
| HCT116 | ensemble | 0.999786 | 0.040293 | 0.029304 |
| HCT116 | specialist | 0.99798 | 0.428571 | 0.355311 |
| HCT116 | student | 0.999157 | 0.0659341 | 0.0695971 |
| HEK293T | ensemble | 0.999252 | 0.0390071 | 0.0851064 |
| HEK293T | specialist | 0.993369 | 0.882979 | 0.939716 |
| HEK293T | student | 0.997123 | 0.0992908 | 0.152482 |
| K562 | ensemble | 0.998582 | 0.0873016 | 0.0793651 |
| K562 | specialist | 0.994061 | 1.13889 | 1.11508 |
| K562 | student | 0.994797 | 0.210317 | 0.206349 |

All response-screen metrics are public summary proxies. Whole-family exclusion prevents using held-family intervention labels in a predictor; previously inspected source biology still limits inferential novelty. Source-specialist weights of zero mean that expert was rejected, not secretly included as an independent vote.

## Disjoint NTC diagnostics

| candidate | predicted_calls | reference_calls | total_CPM_variance | common_expression_numerator |
| --- | --- | --- | --- | --- |
| compact | 156 | 0 | 3.43063e+07 | -21.1981 |
| quantile10 | 235.333 | 0 | 3.16281e+07 | -20.0155 |
| quantile25 | 301 | 0 | 3.49668e+07 | -19.3162 |

Predicted calls compare generated controls with independent right-hand controls. The actual NTC diagnostic cells are disjoint from templates but are included in the left query-mean fitting pool; their expression numerator is optimistic and is not the null selection threshold. A distribution finalist may not increase generated-versus-independent-control DE calls by >10% relative to compact emission. NTC diagnostics use 200 cells per arm, not a 400-cell challenge treatment test.

## Actual-count development

Higher is better for PDS/fidelity/reach/Jaccard; lower for NMAE/expression ratio. Family-balanced NMAE averages only matched eligible panels. Large expression ratios can be dominated by tiny control-distance denominators; the predeclared guard instead compares paired absolute numerators, equivalent to relative ratio worsening when the shared denominator resolves. Missing primary metrics are disclosed rather than filled with zero.

| candidate | pds | fidelity | reach | jaccard | nmae | common_expression_ratio |
| --- | --- | --- | --- | --- | --- | --- |
| response | 0.777778 | 0.604757 | 0.475549 | 0.0221776 | 0.848839 | 68.0079 |
| response_quantile10 | 0.777778 | 0.596487 | 0.526688 | 0.020847 | 0.859735 | 68.0394 |
| response_quantile25 | 0.777778 | 0.588099 | 0.439432 | 0.0200071 | 0.847067 | 68.0661 |
| shuffled_target25 | 0.777778 | 0.593589 | 0.512162 | 0.0210073 | 0.869594 | 68.0439 |
| v13 | 0.777778 | 0.606027 | 0.479518 | 0.0222097 | 0.851085 | 67.7562 |
| v13_quantile25 | 0.777778 | 0.589404 | 0.460723 | 0.0196835 | 0.858126 | 67.7791 |

Available-metric guard comparisons:

| candidate | passes_available_guards | complete_primary | pds | fidelity | reach | jaccard | nmae | relative_expression_worsening |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| response | False | False | 0 | -0.00126953 | -0.00396902 | -3.2036e-05 | -0.00224501 | 0.0223421 |
| v13_quantile25 | False | False | 0 | -0.016623 | -0.0187954 | -0.00252614 | 0.00704138 | 0.029942 |
| response_quantile10 | False | False | 0 | -0.00953967 | 0.0471702 | -0.00136268 | 0.00865023 | 0.0556765 |
| response_quantile25 | False | False | 0 | -0.0179272 | -0.0400861 | -0.00220259 | -0.00401718 | 0.0843209 |

Frozen finalist: **v13**. Equal metric-improvement ranks select among passing candidates; these ranks are not private competition anchor scaling. Full check details and missing-NMAE panels are in `selection.json`. The 2% expression tolerance is an exploratory prospective bound, not statistical noninferiority and not a retroactive change to V14's outcome.

## Matched negative control and reserved verification

Matched distribution control: **not reached**.

The original shuffled25 development arm matches response/quantile25 exactly, but differs in blend from response/quantile10 and in means from V13/quantile25. Before generated-count outcomes, an extra freeze required an exactly mean/blend-matched shuffled-target check for any distribution finalist. It can reject the finalist without retuning or selecting another from the same outcome.

Reserved algorithm verification: **not passed / not reached**.

Reserved verification was not opened after the earlier rejection.

A pass on available guards is not complete six-metric verification when NMAE panels are missing. Current target-name splits and paired seeds do not create new biological lineages. No candidate is reselected from a failed reserved result.

## Integrity and reproduction

Nine scientific tests passed (five new V15 tests plus four imported V14 invariants). The independent experiment audit verifies original unique IDs and disjoint splits, manually reconstructs **2160** source quantile coordinates, checks nested grid choices/cache hashes, and exactly replays **18** original-baseline count/profile/DE cases on V15's exact inputs. The initial historical-hash assertion failed: V14 sampled unsorted control pools before saving sorted pools, whereas V15 sampled the saved pools. Actual treatment identities and whole control halves match, but template identities differ. The failed auditor is preserved as `audit_v15_experiment_initial.py`; `audit_correction.json`, `historical_baseline_differences.csv` and `same_input_baseline_replay.json` record the correction. Frozen methods, V15 identities, outcomes and candidate selection were not changed. This separates numerical integrity from a performance claim. If a package is present, its independent release audit also checks every payload count/index/pointer/axis, every sampled original control identity/guide/library, all production response/distribution arrays and twelve complete-group regenerations.

Root scripts: `v15_ensemble_methods.py`, `v15_population_methods.py`, `prepare_v15_experiment.py`, `run_v15_experiment.py`, `assess_v15_experiment.py`, `check_v15_matched_control.py`, `audit_v15_experiment.py`, `replay_v15_baseline.py`, `audit_v15_experiment_initial.py` (preserved failed audit), `fit_v15_submission.py`, `build_v15_submission.py`, `package_v15_submission.py`, `validate_v15_submission.py`, `release_v15_submission.py`, `test_v15_ensemble.py` and this report generator. Exact commands are in [V15_STRATEGY.md](V15_STRATEGY.md). The release scheduler only sequences frozen stages and stops on rejection; it performs no upload.

Plans, source/model sidecars, scores, selections and audits are under `experiments/v15_weighted_transport/analysis/`. Large arrays, generated counts, packages and full sampling provenance are regenerated locally and excluded from Git. `report_manifest.json` pins the evidence and generated report bytes.
