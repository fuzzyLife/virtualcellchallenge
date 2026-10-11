# V13 submission and validation

Published as `AvAsYkuxFSjbxYDCd0wJ`, overall **0.137883134**, a gain of **+0.000879160** (0.64%) over V11 on the same panel and anchors. This is the best recorded project validation score. Rank was **558** in the user report and **559** at the API check; ranks change with other entries. See [V13_LEADERBOARD_RESULTS.md](V13_LEADERBOARD_RESULTS.md).

**Generation-time assessment:** local count confirmation was inconclusive because some primary NMAE panels were ineligible. The model was packaged as exploratory. This remains the historical validation outcome despite the later small leaderboard gain.

A new package is ready at `experiments/v13_response_programs/prediction.vcc`. Official prep and independent full-payload/replay checks passed. It is now published; the original generation-time report is archived in `experiments/v13_response_programs/analysis/generation_submission_document.md`.

- Selected model: `unknown_kd` from `v13_quality_methods`.
- Package size: 4,081,111,040 bytes (3.80 GiB).
- SHA-256: `2454a8b8ce694217962b89f1bac4ec0fb1bbd01a396c1ffe5cd780f4baa11591`.
- Cells/genes: 360,000 by 18,533; 900 context-target groups with 400 cells each.
- Nonzero counts: 2,375,909,107.
- Prior published reference: V11, overall 0.137003974. The V13 publication and metric comparison are recorded in [V13_LEADERBOARD_RESULTS.md](V13_LEADERBOARD_RESULTS.md).

## Research and tests

See [V13_STRATEGY.md](V13_STRATEGY.md) for the primary literature, complete response alternatives, data coverage, validation design, correction record, and reproduction. We tested six low-rank variants, three robust-consensus variants, and three assay-aware knockdown variants at summary level. Two alternatives from each of the three classes advanced independently to generated-count development. One candidate was frozen before reserved confirmation, then compared once with V11. Development comprises 58 target-seed cases from three source families; confirmation comprises 30 cases from HCT116 and HEK293T only, with no H1 confirmation. Reused source cells, technical strata and paired seeds are not independent biological replication.

### Development (full independent controls)

| candidate | pds | fidelity | reach | jaccard | nmae | common_expression_ratio |
| --- | --- | --- | --- | --- | --- | --- |
| unknown_kd | 0.720833 | 0.586515 | 0.356942 | 0.059946 | 0.878184 | 1.402759 |
| v11 | 0.720833 | 0.586792 | 0.356918 | 0.060104 | 0.877241 | 1.404989 |

### Confirmation (full independent controls)

| candidate | pds | fidelity | reach | jaccard | nmae | common_expression_ratio |
| --- | --- | --- | --- | --- | --- | --- |
| unknown_kd | 0.641667 | 0.618991 | 0.864925 | 0.003983 | 0.784728 | 5.474330 |
| v11 | 0.641667 | 0.617958 | 0.864925 | 0.003972 | 0.784728 | 5.509551 |

Lower NMAE/expression error is better; higher PDS/fidelity/reach/Jaccard is better. These are **unanchored public-source metrics**, not competition scores. Both reference modes are retained in CSV. Some sampled400 panels lack an eligible NMAE contrast; they are reported as missing, never zero. Full independent development panels are complete. Missing confirmation panels, if any, are listed in the assessment records and prevent a full primary-guard pass. Any displayed NMAE average then covers eligible panels only; it is not evidence that the complete confirmation passed. A reporting-only adapter preserves these missing values without altering frozen predictions, metrics, or selection.

Guard results:

- development / sampled400: all guards pass
- development / full_independent: all guards pass
- confirmation / sampled400: all guards pass
- confirmation / full_independent: primary_complete_or_secondary_missing_reported

Conditional target-bootstrap intervals (2,000 Bayesian draws; candidate minus V11):

| stage | metric | delta | low | high | families |
| --- | --- | --- | --- | --- | --- |
| development | fidelity | -0.000277 | -0.000864 | 0.000127 | 3 |
| development | reach | 0.000023 | 0.000001 | 0.000074 | 3 |
| development | jaccard | -0.000158 | -0.000513 | -0.000002 | 3 |
| development | nmae | 0.000942 | 0.000915 | 0.000959 | 3 |
| confirmation | fidelity | 0.001033 | 0.000196 | 0.002192 | 2 |
| confirmation | reach | 0.000000 | 0.000000 | 0.000000 | 2 |
| confirmation | jaccard | 0.000011 | -0.000177 | 0.000195 | 2 |
| confirmation | nmae | 0.000000 | 0.000000 | 0.000000 | 2 |

These intervals condition on the available public source families and do not measure uncertainty over unseen ABC lineages. Missing NMAE panels are excluded explicitly; zero change on an eligible subset is not complete validation. PDS and expression error are nonlinear panel metrics, so their paired family ranges are recorded separately in `selected_uncertainty.json`.

The failed low-rank stage remains recorded and was not relabeled as an improvement. Public development anchored overall: 0.013240 versus V11 0.013840 (difference -0.000600). The raw-guard policy selects expression-error improvement while limiting other losses; this diagnostic does not demonstrate a better overall score.

Official public-source anchored comparisons are separate diagnostics; their small denominators and different contexts prevent interpreting them as leaderboard forecasts. Adaptive development selection increases optimism; frozen count confirmation is the final local check.

## Integrity and reproduction

The held-family mutation tests passed for all three response classes. All V11 baseline profiles, jackknifes and DE metrics replay exactly against independently generated V12 caches. The production response array was independently refitted and compared exactly. A separate post-freeze audit records small float64/float32 pooling differences outside the 65 KD-corrected targets (maximum 8.345e-7 log2 units across 235 other targets); those differences are part of the tested candidate, so this is not a bitwise single-change ablation. The payload audit compares every count, CSR index/pointer and axis; verifies all 360,000 sampled control identities, guides and depths; and exactly regenerates 12 complete groups (88,958,400 gene counts). Counts are nonnegative raw integers, with the required genes/order and cell/group sizes. File integrity does not establish biological prediction accuracy.

After completed development and selection:

```cmd
.\.competition_env\Scripts\python.exe audit_v13_baseline_replay.py
.\.competition_env\Scripts\python.exe select_v13_submission.py
.\.competition_env\Scripts\python.exe evaluate_v13_confirmation.py --workers 4
.\.competition_env\Scripts\python.exe confirm_v13_submission.py
.\.competition_env\Scripts\python.exe fit_v13_submission.py
.\.competition_env\Scripts\python.exe build_v13_submission.py --preflight --workers 4
.\.competition_env\Scripts\python.exe build_v13_submission.py --workers 4
.\.competition_env\Scripts\python.exe validate_v13_submission.py
.\.competition_env\Scripts\python.exe summarize_v13_evidence.py
.\.competition_env\Scripts\python.exe document_v13_submission.py
```

Run in a fresh checkout with the documented source data, or preserve the already frozen plan/model/package rather than overwriting it. Large source/count caches and packages are intentionally excluded from Git; root code and JSON/CSV evidence are tracked.

## Historical upload command (already published)

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v13_response_programs\prediction.vcc -m "Assay-aware knockdown transport v13" --wait
```

Inspect the existing entry with `vcc status AvAsYkuxFSjbxYDCd0wJ`. The historical submit command creates a duplicate entry. V11 remains preserved for comparison.
