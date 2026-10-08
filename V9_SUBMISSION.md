# V9: activation-preserving profile calibration

<!-- V9 publication begin -->

Published as `QBpef7kQNiFa4IFSxNi5`, overall **-0.070447281**, rank 1084 in the user report. V4 remains best at **-0.050673910**. Five scaled metrics improved versus V4, but the fidelity loss outweighed their gains. See [V9_LEADERBOARD_RESULTS.md](V9_LEADERBOARD_RESULTS.md).

<!-- V9 publication end -->

Verified published package: `experiments/v9_activation_preserving/prediction.vcc` (3,295,631,360 bytes).
**360,000 cells x 18,533 genes; 300 targets x 400 cells x three contexts A/B/C.**
Official VCC prep and independent full-payload verification passed. Published leaderboard scores are recorded above; the public-development tables below retain their original meaning.
SHA-256: `b16504d89bf238d9b6fa3c38ba65a8108adb43eab11cbd9c556e853ecf15f8cf`.

## Decision and evidence

Frozen generator: **legacy_fit50**; biological effects: **exact unchanged V4/V7 response arrays**.
Calibration chose the best six-metric public proxy among candidates passing the declared neutral, fidelity, PDS and expression guards. Calibration froze before confirmation. Confirmation improvement: **True**, guard passed **True**, proxy utility difference **+0.019533**. No confirmation-based reselection occurred. These are unanchored public metrics, not a predicted leaderboard overall.

## Why this experiment

The prior validation audit found target-specific control origins inflated constant-profile PDS above its correct 0.5 null. Expression corrections also needed whole-gene scalar jackknifes and signed denominators. V8 fixed those calculations, used control-defined technical strata, target-independent templates, source-lineage exclusion and globally disjoint target calibration/confirmation. Its safe full moment decoder improved some expression metrics but lost supported-target fidelity and did not confirm a proxy improvement. A half-variance decoder produced 968-2,403 false DE calls on unperturbed controls and was rejected. See `V8_SUBMISSION.md` and `TRANSFER_VALIDATION_AUDIT.md`.

V9 retains these validation repairs and tests a different count initialization. First simulate V4's binomial downregulation and zero-aware Poisson activation. Then calibrate mean cell composition and pooled expression by 25%, 50% or 100% toward the full control-population response. Calibration starts from the simulated perturbed cells, preserving activation of genes that were zero in the original template. It preserves each simulated cell's library exactly, rather than imposing untreated depths. Integer rounding and bulk repair are upstream reviewed implementations, with verified realized moment errors. Legacy generation remains an explicit baseline.

The correction strength is selected on new count-validation target genes, excluding every V8 tested target globally and V7 count-tested targets in each held source. All targets were already public source-summary data; these are new generated-cell tests, not pristine blinded biology. HCT116/HEK293T eligibility requires >=150 actual cells in the control-defined stratum; H1 requires >=60. Selection uses metadata and own control expression >5 CPM, never treated expression values. Fixed 400-cell template and reference groups are shared across all targets per context/seed, from disjoint untreated populations. No treated batch proportions enter model inputs.

H1 evaluation excludes direct H1 components and its learned prior; HCT/HEK are absent from V4's response. K562/CD4 supply the transferred fallback. A/B/C still have no measured treated training labels. Consequently the central cross-lineage causal transfer uncertainty remains unresolved; control gene reconstruction neural networks cannot supply those missing intervention responses.

## Calibration: all metrics

| candidate | pds | fidelity | reach | jaccard | nmae | common_expression_ratio | utility |
| --- | --- | --- | --- | --- | --- | --- | --- |
| control_null | 0.500000 | 0.059028 | 0.344246 | 0.166667 | 0.930918 | 4.307320 | -0.208431 |
| legacy | 0.513889 | 0.428023 | 0.575636 | 0.086613 | 0.799665 | 4.368228 | -0.096604 |
| legacy_fit100 | 0.500000 | 0.364942 | 0.599069 | 0.087904 | 0.816072 | 2.678414 | -0.035185 |
| legacy_fit25 | 0.513889 | 0.448221 | 0.587733 | 0.065705 | 0.802907 | 3.613268 | -0.083694 |
| legacy_fit50 | 0.500000 | 0.433340 | 0.593906 | 0.065643 | 0.807450 | 3.076719 | -0.060800 |

## Confirmation: all metrics

| candidate | pds | fidelity | reach | jaccard | nmae | common_expression_ratio | utility |
| --- | --- | --- | --- | --- | --- | --- | --- |
| control_null | 0.500000 | 0.083333 | 0.553149 | 0.416667 | 0.927965 | 4.106230 | -0.152049 |
| legacy | 0.493056 | 0.479674 | 0.701401 | 0.094193 | 0.732988 | 4.247059 | -0.068406 |
| legacy_fit100 | 0.500000 | 0.404264 | 0.763848 | 0.182774 | 0.746220 | 3.134337 | -0.029550 |
| legacy_fit25 | 0.500000 | 0.440863 | 0.705285 | 0.097936 | 0.735158 | 3.703373 | -0.061127 |
| legacy_fit50 | 0.506944 | 0.397573 | 0.700885 | 0.143570 | 0.739579 | 3.334415 | -0.048873 |

Raw fidelity follows the official direction-yield budget, not sign precision. Undefined zero-budget values remain undefined. The safety guard uses the same reference-call>=10 targets for each candidate and count-weighted fidelity within 0.02 of legacy; expression <=1.02 legacy, PDS within 0.05, and every neutral-control false-DE fraction <=0.005. This stability guard is a development criterion, not a replacement official metric. Utility is the V8 equal-family average of centered PDS, fidelity, reach, Jaccard, clipped 1-NMAE and clipped 1-expression. Raw expression denominators can be weak, and absent hidden anchor sets prevent translating these values into leaderboard scores. Reference-supported fidelity and guards are recorded separately in the JSON decisions.

## Reference-supported fidelity guard

| stage | candidate | supported_budget_fidelity |
| --- | --- | --- |
| calibration | control_null | 0.000000 |
| calibration | legacy | 0.235849 |
| calibration | legacy_fit100 | 0.215135 |
| calibration | legacy_fit25 | 0.227401 |
| calibration | legacy_fit50 | 0.220419 |
| confirmation | control_null | 0.000000 |
| confirmation | legacy | 0.002476 |
| confirmation | legacy_fit100 | 0.002724 |
| confirmation | legacy_fit25 | 0.002229 |
| confirmation | legacy_fit50 | 0.002724 |

Passing a tolerance guard is not strong absolute causal accuracy. Confirmation supported-target fidelity is very low for both legacy and calibrated transfer; this remains a major response-model weakness. The improved decoder does not supply missing matched A/B/C intervention data or establish that the V4 response directions are correct. The direct target retains V4's fixed 90% knockdown assumption; it is excluded from scored gene sets but still affects library normalization. This experiment does not learn hidden dose or timing. Raw mean fidelity can be dominated by tiny reference budgets and differs substantially from this supported-target summary.

## Metric availability

| stage | family | candidate | panels | finite_utility_panels | resolved_expression_panels | finite_reach_panels |
| --- | --- | --- | --- | --- | --- | --- |
| calibration | H1 | control_null | 2 | 2 | 2 | 2 |
| calibration | H1 | legacy | 2 | 2 | 2 | 2 |
| calibration | H1 | legacy_fit100 | 2 | 2 | 2 | 2 |
| calibration | H1 | legacy_fit25 | 2 | 2 | 2 | 2 |
| calibration | H1 | legacy_fit50 | 2 | 2 | 2 | 2 |
| calibration | HCT116 | control_null | 4 | 4 | 4 | 4 |
| calibration | HCT116 | legacy | 4 | 4 | 4 | 4 |
| calibration | HCT116 | legacy_fit100 | 4 | 4 | 4 | 4 |
| calibration | HCT116 | legacy_fit25 | 4 | 4 | 4 | 4 |
| calibration | HCT116 | legacy_fit50 | 4 | 4 | 4 | 4 |
| calibration | HEK293T | control_null | 4 | 4 | 4 | 4 |
| calibration | HEK293T | legacy | 4 | 4 | 4 | 4 |
| calibration | HEK293T | legacy_fit100 | 4 | 4 | 4 | 4 |
| calibration | HEK293T | legacy_fit25 | 4 | 4 | 4 | 4 |
| calibration | HEK293T | legacy_fit50 | 4 | 4 | 4 | 4 |
| confirmation | H1 | control_null | 2 | 2 | 2 | 2 |
| confirmation | H1 | legacy | 2 | 2 | 2 | 2 |
| confirmation | H1 | legacy_fit100 | 2 | 2 | 2 | 2 |
| confirmation | H1 | legacy_fit25 | 2 | 2 | 2 | 2 |
| confirmation | H1 | legacy_fit50 | 2 | 2 | 2 | 2 |
| confirmation | HCT116 | control_null | 4 | 2 | 2 | 2 |
| confirmation | HCT116 | legacy | 4 | 2 | 2 | 2 |
| confirmation | HCT116 | legacy_fit100 | 4 | 2 | 2 | 2 |
| confirmation | HCT116 | legacy_fit25 | 4 | 2 | 2 | 2 |
| confirmation | HCT116 | legacy_fit50 | 4 | 2 | 2 | 2 |
| confirmation | HEK293T | control_null | 4 | 3 | 3 | 4 |
| confirmation | HEK293T | legacy | 4 | 3 | 3 | 4 |
| confirmation | HEK293T | legacy_fit100 | 4 | 3 | 3 | 4 |
| confirmation | HEK293T | legacy_fit25 | 4 | 3 | 3 | 4 |
| confirmation | HEK293T | legacy_fit50 | 4 | 3 | 3 | 4 |

Undefined metrics stay undefined: reach requires a nonempty real DE budget, and the expression ratio requires a resolved positive corrected denominator. Pandas aggregation averages defined panels/strata/families separately for each metric. Thus a nominal equal-family proxy may include fewer than three families when an entire family's utility is undefined. The table and raw panel CSVs expose those omissions; summary means must not be interpreted as complete three-lineage validation or official anchors. Calibration and confirmation are different targets and cannot be compared as the same test.

## Neutral checks

| generator | max_false_de_calls | mean_false_de_calls | max_false_call_fraction | mean_variance_ratio | neutral_guard |
| --- | --- | --- | --- | --- | --- |
| legacy | 1 | 0.166667 | 0.000099 | 0.987813 | True |
| legacy_fit100 | 1 | 0.166667 | 0.000104 | 1.017966 | True |
| legacy_fit25 | 1 | 0.166667 | 0.000099 | 0.992696 | True |
| legacy_fit50 | 1 | 0.166667 | 0.000104 | 1.000327 | True |

## Coverage

| context | left_controls | right_controls | calibration_targets | confirmation_targets | actual_cells_min | actual_cells_max |
| --- | --- | --- | --- | --- | --- | --- |
| H1_controlstratum0 | 768 | 768 | 2 | 2 | 94 | 200 |
| HCT116_controlstratum0 | 960 | 960 | 4 | 4 | 152 | 200 |
| HCT116_controlstratum1 | 782 | 782 | 4 | 4 | 150 | 192 |
| HEK293T_controlstratum0 | 2144 | 2144 | 4 | 4 | 157 | 200 |
| HEK293T_controlstratum1 | 1414 | 1414 | 4 | 4 | 153 | 200 |

Public raw release selection contains sampled HCT/HEK treated cells, 25 H1 source targets and incomplete K562 gene coverage. Existing source provenance records quantify these limitations. Technical strata control batches; they do not identify hidden context lineage, timing, knockdown efficiency or platform. A confirmed public proxy improvement would still require a leaderboard evaluation. A failed confirmation is preserved openly and cannot be interpreted as demonstrated improvement.

## Reproduce

Use the existing pinned source caches, reviewed upstream clone and completed V8 frozen plan/model. `run_v8_pipeline.py --resume` reconstructs its prerequisites before V9. Root scripts retain every step. To reproduce from scratch, preserve the frozen experiment directory under a new name first; scripts refuse to overwrite frozen decisions. Then:

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe run_v9_pipeline.py --workers 4
```

`--resume` preserves completed stages. The pipeline runs scientific tests, reserves fresh targets, evaluates calibration, tests disjoint challenge NTCs, freezes selection, evaluates confirmation once, freezes production dependencies, audits splits, preflights 24 groups, generates 900 groups, runs official prep and independently checks the full payload. It does not upload.

The independent audit checks all 360,000 source IDs, guide/depth records and integer count libraries; compares every packaged sparse value/index/pointer and axis against the generated matrix; and replays 12 complete groups (88,958,400 gene counts) exactly. Model and package SHA-256 plus code dependency hashes are saved in `analysis/`. `audit_v9_provenance.py` separately verifies linked calibration/confirmation/neutral record hashes, source file hashes, prior model references and all V7/V8 target exclusions without changing frozen numerical inputs. `.gitattributes` preserves the exact code and experiment-record bytes rather than rewriting line endings.

## Historical upload command (already published)

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v9_activation_preserving\prediction.vcc -m "Activation-preserving profile calibration v9" --wait
```

This upload is complete. Inspect the existing result with `vcc status QBpef7kQNiFa4IFSxNi5`. Repeating the submit command creates a duplicate submission.
