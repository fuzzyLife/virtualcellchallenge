# V19 results: covariance corrections remain a research result

Completed 2026-10-13. **The covariance candidate failed the predeclared release criteria. No V19 submission was built or uploaded. Reserved verification remains unopened.** V16 remains the released reference.

The experiment used five arms and the same 28 paired development cases as V18: **140 full-count population comparisons**. Each uses 400 distinct original template cells, 400 genuine treatment cells and the same independent right controls. All 300 panel genes are excluded from the training program basis and PDS; the own target is excluded from DE scoring. V16 response fitting and the new program/force fitting exclude the scored source family completely.

| Raw metric | V16 | Covariance | Diagonal control | Amplitude control | Better |
|---|---:|---:|---:|---:|---|
| PDS | 0.777778 | 0.794444 | 0.766667 | 0.777778 | Higher |
| Fidelity | 0.606476 | 0.607318 | 0.606121 | 0.607092 | Higher |
| Reach | 0.503430 | 0.501982 | 0.427019 | 0.506239 | Higher |
| Jaccard | 0.022365 | 0.022333 | 0.022003 | 0.022682 | Higher |
| Eligible-panel NMAE | 0.849526 | 0.853512 | 0.856864 | 0.854622 | Lower |

Covariance improved PDS by 0.01667 and fidelity by 0.00084. Reach/Jaccard declined slightly; NMAE worsened by 0.00399 and paired expression numerator by 0.296%. Those declines stay inside the practical guards, but only two metrics improved, versus the required three. The model beat diagonal and shuffled covariance on two DE measures, but failed to beat the amplitude-matched comparator on two DE measures. The apparent benefit therefore does not justify replacing V16. These are raw public metrics, not private anchor-scaled leaderboard scores.

HCT116 NMAE remains unavailable in both seeds because true DE counts are insufficient. The public development biology has been used before; repeated draws and source-family holds do not create new independent biological lineages. No unseen-private-context performance or ranking gain is claimed.

## Independent checks

Mathematical unit tests passed: identity/scalar covariance maps; singular-value bounds; nuisance-depth centering; zero-force fallback; correction bounds; measured-gene masking; low-coverage abstention; own-target preservation.

A separate dense-arithmetic audit independently reproduced **nine source control covariance matrices** across the three excluded-family models and all-source model. It also reproduced eligible source treatment shifts for fixed sampled targets (ABCD1, HDAC8, VIM, SMARCA5). It used direct least squares and cross-products rather than the fitted sparse pipeline. Source cache SHA-256, actual source cell counts, basis orthogonality, positive covariance eigenvalues, panel exclusion, case identities and output hashes passed. All 28 V16 comparator count digests again matched the independently tested V17 comparator.

Source support: H1 has 25 panel targets with at least 40 selected cells; HCT116 275; HEK293T 287. The held-family basis uses only the other sources' control variability. Query corrections use only frozen NTC rows; no query treatment amplitude or response is supplied. See [V19_STRATEGY.md](V19_STRATEGY.md) for exact equations, bounds and reproduction commands.

Evidence is under `experiments/v19_covariance_forces/analysis/`: `plan.json`, `models/*.json`, `model_input_independent_audit.json`, `development_independent_audit.json`, `development_summary.csv`, `development_panels.csv`, `development_scores.csv`, `development_expert_diagnostics.csv`, and `selection.json`.

The next approved research path is a bounded paired-control Stack population-model pilot. It will not inherit a release claim from these small covariance gains.
