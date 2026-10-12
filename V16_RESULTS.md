# V16 measured results: control-preserving confidence transport

V16 is published as `wG3GGbaTsl2dXFn7cQy1`, overall **0.141613857**, a gain of **+0.003730723** (2.71%) over V13 on the same panel and anchors. It is the best recorded project validation score. Rank was **562** in the user report and **561** at the API check; historical ranks are not directly comparable. See [V16_LEADERBOARD_RESULTS.md](V16_LEADERBOARD_RESULTS.md).

A new exploratory package passed official prep and full-payload validation. At generation time it had not been uploaded and no score was known; publication is recorded below.

V13 was the best published project score at generation time (0.13788313368012328). Public raw comparisons below are not an official anchored score or a prediction of rank. Current public source biology and development targets are partly reused; missing NMAE panels prevent complete six-metric verification.

## Implemented and tested

Three response forms (V13, gene-calibrated shrinkage, matched uniform shrinkage) were crossed with compact counts, balanced molecule redistribution and residual10. The fixed 18-bin calibrator uses inner-family held predictions, excludes the outer family and all reserved calibration-label targets, and retains at least 50% V13. Two-source significant effects >=0.2 log2 remain intact. CD4 conditions are one family; source significance is an approximate confidence feature, not significance in private A/B/C contexts.

Balanced transport changes each sampled original cell by binomial removals/multinomial redistribution, preserves libraries exactly and is bitwise identity under a zero intervention. It preserves original zero support and cannot activate a gene absent in an individual cell. It does not force global population moments. Residual10 retains the original 0.65 Poisson coefficient before fitting; the dual solver can still change final dispersion. These are exploratory generators, not optimal transport or a causal biological model. Full equations, identities, thresholds and commands are in [V16_STRATEGY.md](V16_STRATEGY.md).

Nine scientific unit tests passed. Counts/comparisons: **{'null': {'population_comparisons': 18, 'newly_generated_populations': 18, 'reused_populations': 0, 'target_seed_cases': 6}, 'development': {'population_comparisons': 252, 'newly_generated_populations': 224, 'reused_populations': 28, 'target_seed_cases': 28}, 'confirmation': {'population_comparisons': 36, 'newly_generated_populations': 36, 'reused_populations': 0, 'target_seed_cases': 18}}**. These counts include repeated seeds, shared truths and reused baselines; they are not independent biological experiments. The actual development treatments contain 400 distinct cells each. Null comparisons generate 400 controls against independent right controls; their disjoint actual-control diagnostic contains 200 cells. Source validation shares templates within a family/seed; production samples controls separately per target. No hidden competition response labels are used.

## Response calibration diagnostics

These are summary proxies and did not choose candidates; actual-count gates determine release.

| candidate | mse_ratio | fidelity_proxy | jaccard_proxy | selection_objective |
| --- | --- | --- | --- | --- |
| calibrated | 1.02357 | 0.212997 | 0.0838801 | 0.121575 |
| v13 | 1.0798 | 0.278998 | 0.0664159 | 0.136248 |

| excluded_family | minimum_gain | median_gain | fully_retained_bins | training_gene_observations |
| --- | --- | --- | --- | --- |
| none (production) | 0.0847667 | 0.938762 | 7 | 14429632 |
| CD4 | 0.0965651 | 0.938119 | 3 | 8660596 |
| H1 | 0.0797372 | 0.967682 | 8 | 14156814 |
| HCT116 | 0.0843695 | 0.859388 | 5 | 11345541 |
| HEK293T | 0.0641555 | 0.918934 | 7 | 11003345 |
| K562 | 0.0842781 | 0.934325 | 6 | 12552232 |

Gene-observation counts contain correlated genes, targets and conditions; they are not statistical sample sizes. Gains are regularized toward one and clipped, with explicit strong-evidence protection. Uniform controls match each calibrated target's total absolute log2 effect; they remain eligible candidates under the same prospective policy.

## Independent-control diagnostics

| candidate | predicted_calls | reference_calls | variance_ratio | variance_log_distance | detection_MAE |
| --- | --- | --- | --- | --- | --- |
| balanced | 0 | 0 | 0.975819 | 0.0276389 | 0.0166775 |
| compact | 1699.17 | 0 | 0.196021 | 1.63614 | 0.0716696 |
| residual10 | 1522.83 | 0 | 0.205201 | 1.58989 | 0.0699183 |

Zero false calls in a control-only diagnostic establish neither perturbation accuracy nor private performance. The generator gate requires <=10% more false calls than compact, improved absolute log total-variance distance, and detection MAE within 10% of compact. Total variance alone does not verify every gene or covariance.

## Paired actual-count development

Higher is better for PDS/fidelity/reach/Jaccard; lower is better for NMAE/expression. The displayed expression ratio can be dominated by tiny control-distance denominators; the policy uses paired expression numerators for relative worsening. NMAE averages eligible matched panels; missing panels are disclosed in `selection.json`, never filled with zero.

| candidate | pds | fidelity | reach | jaccard | nmae | common_expression_ratio |
| --- | --- | --- | --- | --- | --- | --- |
| balanced | 0.752778 | 0.283474 | 0.496519 | 0.169599 | 0.862612 | 257.883 |
| calibrated_balanced | 0.719444 | 0.207576 | 0.491473 | 0.170336 | 0.804239 | 256.767 |
| calibrated_compact | 0.713889 | 0.608289 | 0.424309 | 0.0216081 | 0.80514 | 66.9599 |
| calibrated_residual10 | 0.713889 | 0.610509 | 0.446047 | 0.0220349 | 0.793717 | 66.9449 |
| residual10 | 0.777778 | 0.606476 | 0.50343 | 0.0223651 | 0.849526 | 67.724 |
| uniform_balanced | 0.719444 | 0.194552 | 0.488794 | 0.172948 | 0.794479 | 256.953 |
| uniform_compact | 0.713889 | 0.607323 | 0.500085 | 0.0212385 | 0.782696 | 66.9631 |
| uniform_residual10 | 0.713889 | 0.609454 | 0.492482 | 0.0216565 | 0.793647 | 66.9373 |
| v13 | 0.777778 | 0.606027 | 0.479518 | 0.0222097 | 0.851085 | 67.7562 |

| candidate | passes | pds | fidelity | reach | jaccard | nmae | relative_expression_worsening |
| --- | --- | --- | --- | --- | --- | --- | --- |
| calibrated_compact | False | -0.0638889 | 0.00226255 | -0.0552095 | -0.00060158 | -0.0459442 | -0.0417947 |
| balanced | False | -0.025 | -0.322553 | 0.0170003 | 0.14739 | 0.0115275 | 1.76121 |
| calibrated_balanced | False | -0.0583333 | -0.398451 | 0.0119549 | 0.148127 | -0.0468457 | 1.754 |
| residual10 | True | 0 | 0.000449189 | 0.0239115 | 0.000155425 | -0.00155884 | -0.000916606 |
| calibrated_residual10 | False | -0.0638889 | 0.00448219 | -0.0334708 | -0.000174797 | -0.0573675 | -0.0418775 |
| uniform_compact | False | -0.0638889 | 0.00129694 | 0.0205668 | -0.000971171 | -0.0683889 | -0.026432 |
| uniform_balanced | False | -0.0583333 | -0.411474 | 0.00927561 | 0.150738 | -0.0566058 | 1.76621 |
| uniform_residual10 | False | -0.0638889 | 0.00342722 | 0.0129638 | -0.000553144 | -0.0574376 | -0.0265844 |

Frozen finalist: **residual10**. If the finalist is `residual10`, deployment keeps the exact V13 response predictor and changes only count generation; the gene-calibrated shrinkage is not deployed. Every candidate must meet the unchanged V15 raw/worst-family bounds and improve at least three available metrics including a DE metric. Gene-calibrated arms must additionally beat their uniform counterparts on two of fidelity/reach/Jaccard. Equal improvement ranks select among passing candidates; no private anchor utility is estimated.

## Reserved verification and release

| candidate | pds | fidelity | reach | jaccard | nmae | common_expression_ratio |
| --- | --- | --- | --- | --- | --- | --- |
| residual10 | 0.611111 | 0.605079 | 0.523197 | 0.0228178 | 0.94183 | 11.2949 |
| v13 | 0.611111 | 0.608499 | 0.480207 | 0.0232707 | 0.94487 | 11.3032 |

Paired reserved changes (candidate minus V13; expression is relative numerator worsening):

| metric | change |
| --- | --- |
| pds | 0 |
| fidelity | -0.00341985 |
| reach | 0.0429894 |
| jaccard | -0.000452901 |
| nmae | -0.00304047 |
| relative_expression_worsening | -0.000633226 |

Reserved result: **True**. This is one frozen finalist on nine reserved current target names, with its uniform control if gene-calibrated. It is not independent new biology. Failure cannot retune or reselect. The plan and once-only assessments are hashed and preserved.

## What the comparison teaches

Balanced transport fixed the control-only artifact but lost 0.323 development fidelity and increased paired expression error by 176.1%. It was rejected despite much better DE overlap. These tests show that preserving control variability alone did not recover perturbation accuracy; sample means and transferred responses still matter.

Gene-calibrated compact responses improved eligible NMAE by 0.0459 and expression error by 4.18%, but reduced PDS by 0.0639 and reach by 0.0552. Shrinking weak evidence also removed useful target separation. It was rejected.

Residual10 reduced null DE calls by 10.4%, but retained only 20.5% of reference total variance. This is a modest count-generation compromise; the major distribution mismatch remains unresolved. The available public checks cannot establish a significant leaderboard improvement.

## Independent integrity audit

The audit verifies unique original IDs, disjoint controls, exact copied actual/template/reference identities, manual gain equations and inner/outer exclusions, all completed count/profile cache hashes, **28** exact reused V15 baseline count/profile/DE cases, **6** baseline regenerations through the original V11 emitter, and **six** exact balanced-control identities. It also verifies the earlier V15 evidence is unchanged. Format/integrity is distinct from biological accuracy.

Final release validation passed: every 2,368,784,027 nonzero count and index, every pointer and axis, all 360,000 original-control identities/guides/libraries, all 900 groups, and 12 exactly regenerated groups (88,958,400 gene-count coordinates). Production responses were recomputed exactly. Package SHA-256: `aec15ac5e5ef38e835bb4206f62c431772a28690fe5823490f045c5f0319cc08`.

The source file and official staging file were compared losslessly before official prep. Balanced production, if selected, additionally checks every predicted nonzero was present in its original cell; the released residual10 generator does not enforce that support restriction.

Root code includes `v16_methods.py`, `test_v16_transport.py`, `prepare_v16_experiment.py`, `run_v16_experiment.py`, `assess_v16_experiment.py`, `audit_v16_experiment.py`, `fit_v16_submission.py`, `build_v16_submission.py`, `package_v16_submission.py`, `validate_v16_submission.py`, `release_v16_submission.py` and this report generator. Small evidence is versioned in `experiments/v16_control_preserving/analysis/`; large arrays/counts/packages/provenance are regenerated locally and excluded from Git. The generation pipeline performs no upload; the user subsequently uploaded this package.
