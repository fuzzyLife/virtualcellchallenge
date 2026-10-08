# V8 submission: shared-control validation and calibrated count profiles

Upload-ready file: `experiments/v8_corrected_transfer/prediction.vcc` (3,275,960,320 bytes). It contains **360,000 cells x 18,533 genes**, all **300 targets x 400 cells x A/B/C**. Official VCC prep and independent full-payload validation passed. This package has not been uploaded and has no leaderboard score.

SHA-256: `fe8d89ca5b147a4d2208a8eac164058a1a892ead45eb1afc1702dd3dbd724dd7`.

## Final decision

Candidate **baseline_dual**, response **v4_a100**, count generator **dual**. The response profile screen favored `v4_a65`; the generated-cell safety selection chose the candidate above. Confirmation improvement under the declared public proxy: **False**, utility difference **-0.003444**. Confirmation was assessed once without reselection. Public proxy utility is not an official anchored leaderboard score.

The original count-only calibration choice `selected_dual_half` was not accepted automatically. An independent neutral-control test rejected `dual_half`: it produced 968-2,403 significant genes without any applied perturbation. Both `legacy` and `dual` produced only 0-1 false calls across the six control tests. V8 uses the safe frozen calibration winner.

The response-effect arrays are exactly equal to the V4/V7 fallback arrays (maximum difference 0). This is an isolated count-decoder change, not a new biological response model. Unlike the V7 resampling comparison, it changes the expected mean/pooled profiles and repairs both integer count moments.

## What was fixed

- One common reference control profile and cell group per evaluation stratum/seed replaces target-specific control origins. Constant predictions now score exactly 0.5.
- Control-only batch PCA/KMeans defines technical strata before response modeling. Prediction templates use a target-independent 400-cell draw, never actual treated batch proportions.
- Source lineages are excluded from their own response training. Calibration and confirmation target genes are globally disjoint. None of these selected targets were used in V7 generated-cell tests; they remain public development data already available to source-summary work.
- Expression calculations retain signed corrections, scalar whole-gene jackknifes, and explicit unresolved-denominator reporting. DE uses the same common control, its CPM>5 gene universe, and per-target BH correction.
- Response pooling, shrinkage and amplitude were screened separately from the count decoder. Full moment matching anchors mean cell composition and pooled expression to the entire available control population, preserves cell depths, and verifies realized integer moments.
- Neutral/control and independent-replicate diagnostics are recorded. A reduced-variance decoder was rejected for creating extensive false DE rather than adopted for its apparent metric gain.

## Public calibration

| candidate | pds | fidelity | reach | jaccard | nmae | common_expression_ratio | utility |
| --- | --- | --- | --- | --- | --- | --- | --- |
| baseline_dual | 0.520833 | 0.205888 | 0.651204 | 0.263944 | 0.881138 | 6.173662 | -0.004656 |
| baseline_legacy | 0.534722 | 0.400271 | 0.656502 | 0.074979 | 0.876142 | 10.925247 | -0.080185 |
| control_null | 0.500000 | 0.045833 | 0.592606 | 0.312500 | 0.915827 | 10.270611 | -0.123003 |
| selected_dual | 0.520833 | 0.054702 | 0.649449 | 0.344055 | 0.893644 | 5.741509 | -0.012541 |
| selected_dual_half | 0.520833 | 0.643933 | 0.475632 | 0.011726 | 0.895677 | 5.720941 | 0.001277 |
| selected_legacy | 0.520833 | 0.274666 | 0.634426 | 0.220730 | 0.891490 | 10.498985 | -0.083181 |

## Public confirmation

| candidate | pds | fidelity | reach | jaccard | nmae | common_expression_ratio | utility |
| --- | --- | --- | --- | --- | --- | --- | --- |
| baseline_dual | 0.611111 | 0.309739 | 0.659669 | 0.131738 | 0.777922 | 5.819636 | 0.087160 |
| baseline_legacy | 0.604167 | 0.534977 | 0.678963 | 0.102908 | 0.772291 | 10.257109 | 0.090604 |
| control_null | 0.500000 | 0.067708 | 0.449298 | 0.208333 | 0.940138 | 9.961570 | -0.072790 |
| selected_dual | 0.569444 | 0.019033 | 0.667450 | 0.238178 | 0.840314 | 5.519274 | 0.034194 |
| selected_dual_half | 0.569444 | 0.669933 | 0.537313 | 0.018744 | 0.837236 | 5.503214 | 0.085038 |
| selected_legacy | 0.569444 | 0.198108 | 0.655157 | 0.222539 | 0.827988 | 9.980427 | 0.030150 |

Every metric is shown, including tradeoffs. Fidelity is the raw yield metric, not sign precision. Many HCT116/HEK293T development targets have zero or tiny reference DE sets, which makes unweighted fidelity comparisons and candidate-specific undefined values unstable. Before confirmation, the final safety selector retained expression/PDS guards, required every neutral false-DE fraction <=0.5%, and replaced its initial mean-fidelity guard with a common reference-supported budget-weighted fidelity guard (reference calls >=10; within 0.02 of baseline). The original calibration decision and all raw metric means remain preserved. This is a calibration refinement, not a claim that the budget-weighted guard is the official leaderboard metric.

## Metric availability

| stage | family | candidate | panels | finite_utility_panels | resolved_expression_panels | finite_reach_panels |
| --- | --- | --- | --- | --- | --- | --- |
| calibration | H1 | baseline_dual | 2 | 2 | 2 | 2 |
| calibration | H1 | baseline_legacy | 2 | 2 | 2 | 2 |
| calibration | H1 | control_null | 2 | 2 | 2 | 2 |
| calibration | H1 | selected_dual | 2 | 2 | 2 | 2 |
| calibration | H1 | selected_dual_half | 2 | 2 | 2 | 2 |
| calibration | H1 | selected_legacy | 2 | 2 | 2 | 2 |
| calibration | HCT116 | baseline_dual | 4 | 4 | 4 | 4 |
| calibration | HCT116 | baseline_legacy | 4 | 4 | 4 | 4 |
| calibration | HCT116 | control_null | 4 | 4 | 4 | 4 |
| calibration | HCT116 | selected_dual | 4 | 4 | 4 | 4 |
| calibration | HCT116 | selected_dual_half | 4 | 4 | 4 | 4 |
| calibration | HCT116 | selected_legacy | 4 | 4 | 4 | 4 |
| calibration | HEK293T | baseline_dual | 4 | 4 | 4 | 4 |
| calibration | HEK293T | baseline_legacy | 4 | 4 | 4 | 4 |
| calibration | HEK293T | control_null | 4 | 4 | 4 | 4 |
| calibration | HEK293T | selected_dual | 4 | 4 | 4 | 4 |
| calibration | HEK293T | selected_dual_half | 4 | 4 | 4 | 4 |
| calibration | HEK293T | selected_legacy | 4 | 4 | 4 | 4 |
| confirmation | H1 | baseline_dual | 2 | 2 | 2 | 2 |
| confirmation | H1 | baseline_legacy | 2 | 2 | 2 | 2 |
| confirmation | H1 | control_null | 2 | 2 | 2 | 2 |
| confirmation | H1 | selected_dual | 2 | 2 | 2 | 2 |
| confirmation | H1 | selected_dual_half | 2 | 2 | 2 | 2 |
| confirmation | H1 | selected_legacy | 2 | 2 | 2 | 2 |
| confirmation | HCT116 | baseline_dual | 4 | 4 | 4 | 4 |
| confirmation | HCT116 | baseline_legacy | 4 | 4 | 4 | 4 |
| confirmation | HCT116 | control_null | 4 | 4 | 4 | 4 |
| confirmation | HCT116 | selected_dual | 4 | 4 | 4 | 4 |
| confirmation | HCT116 | selected_dual_half | 4 | 4 | 4 | 4 |
| confirmation | HCT116 | selected_legacy | 4 | 4 | 4 | 4 |
| confirmation | HEK293T | baseline_dual | 4 | 4 | 4 | 4 |
| confirmation | HEK293T | baseline_legacy | 4 | 4 | 4 | 4 |
| confirmation | HEK293T | control_null | 4 | 4 | 4 | 4 |
| confirmation | HEK293T | selected_dual | 4 | 4 | 4 | 4 |
| confirmation | HEK293T | selected_dual_half | 4 | 4 | 4 | 4 |
| confirmation | HEK293T | selected_legacy | 4 | 4 | 4 | 4 |

Undefined values remain undefined. Aggregation averages defined panels/strata/families separately for each metric; if a family's utility is entirely undefined, its missing value is skipped. A nominal equal-family proxy can therefore include fewer than three lineages. The raw CSVs and availability table expose those omissions. A positive corrected expression denominator and a nonempty real DE budget are required for the corresponding ratio/reach metrics.

## Neutral decoder checks

| generator | max_false_de_calls | mean_false_de_calls | max_false_call_fraction | mean_variance_ratio | neutral_guard |
| --- | --- | --- | --- | --- | --- |
| dual | 1 | 0.166667 | 0.000104 | 1.017966 | True |
| dual_half | 2403 | 1679.500000 | 0.242678 | 0.324835 | False |
| legacy | 1 | 0.166667 | 0.000099 | 0.987813 | True |

## Validation coverage and limits

| context | left_controls | right_controls | calibration_targets | confirmation_targets | observed_cells_min | observed_cells_max | old_v7_count_targets_reused |
| --- | --- | --- | --- | --- | --- | --- | --- |
| H1_controlstratum0 | 768 | 768 | 4 | 4 | 96 | 96 | 0 |
| HCT116_controlstratum0 | 960 | 960 | 4 | 4 | 83 | 145 | 0 |
| HCT116_controlstratum1 | 782 | 782 | 4 | 4 | 62 | 141 | 0 |
| HEK293T_controlstratum0 | 2144 | 2144 | 4 | 4 | 61 | 200 | 0 |
| HEK293T_controlstratum1 | 1414 | 1414 | 4 | 4 | 89 | 187 | 0 |

The source tests contain 60-200 actual cells per target, with 400 generated and reference cells. Technical strata reduce batch variation but cannot establish matching platform, timing, knockdown efficiency or response biology for hidden A/B/C treatments. Replicate diagnostics use independent treated halves and independent control origins; their smaller group sizes are not official anchors. Public scores, including Jaccard gains on empty-response targets, must not be presented as hidden performance estimates. The frozen safe candidate did not satisfy the confirmation improvement rule. This upload package is an experimental calibration candidate, not a demonstrated improvement over V4.

## Reproduce

From Windows Command Prompt, with existing pinned source/control caches:

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe run_v8_pipeline.py --resume --workers 4
```

For a fresh reproduction, omit `--resume` in a checkout without V8 outputs. Stages are: unit tests, control-only plan, response calibration, generated-cell calibration, neutral decoder checks, safety selection, confirmation, final confirmation assessment, production fit, model audit, preflight, full count generation, official prep, independent validation, documentation. Original V4-V7 numerical modules and payloads are preserved. A fresh numerical rerun requires a workspace without the V8 experiment directory; archive the recorded V8 directory before starting rather than overwriting its frozen records.

Key root scripts: `prepare_v8_validation.py`, `v8_methods.py`, `evaluate_v8_submission.py`, `benchmark_v8_neutral.py`, `select_v8_safe_generator.py`, `confirm_v8_safe_selection.py`, `fit_v8_submission.py`, `build_v8_submission.py`, `validate_v8_submission.py`, `test_v8_submission.py`, `run_v8_pipeline.py`, `document_v8_submission.py`.

## Submit

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v8_corrected_transfer\prediction.vcc -m "Shared-control calibrated count profiles v8" --wait
```

If interrupted after upload/scoring starts, use `vcc status <entry-id> --wait` instead of creating another entry.

## Integrity records

`submission_ready.json` holds the package checksum and official prep metadata. `independent_validation.json` confirms every packaged count/index/pointer and axis matches H5AD, scans all 360,000 cell libraries, checks all control IDs and guide sampling, and exactly replays 12 groups (88,958,400 gene-count entries). `reference_model.json` freezes model arrays, code and inputs; selections and evaluation files preserve calibration/confirmation chronology and sampling. `analysis/v8_official_metric_review.json` pins the official source examined for metric parity.
