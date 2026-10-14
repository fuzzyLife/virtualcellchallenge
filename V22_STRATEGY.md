# V22: state-dependent continuous-response pilot

This experiment implements the first bounded pilot in POST_V20_RETROSPECTIVE.md. V16 remains the published champion until a new model earns release. The frozen plan is experiments/v22_state_response/analysis/plan.json. No automatic upload is performed.

## Diagnostic before release

Nine public contrasts cross predicted versus measured treatment-calibration moments with compact versus empirical treated-cell shape, plus a real calibration-cell replica. Calibration and truth are disjoint 200-cell draws: existing caches cap HCT116/HEK293T at 400 cells per target and H1 at about 768, so 400+400 is unsupported. Every diagnostic arm uses calibration treatment depths, which are not available for the competition query. Oracle arms are diagnostic only, never deployable. Mean/bulk fitting errors are measured rather than assumed exact.

## Model

Four control-derived state factors use 768 training-only stable genes; all 300 competition panel genes are excluded. Gene stability is selected using only training-family intervention summaries. Output interactions cover 2,048 measured non-panel genes. Query controls supply their own centering/scaling; source controls are matched by source batch. There is no held-family treatment input during fitting or query inference.

Training-cell residuals are decomposed into a shared across-target component and a target-specific fingerprint after subtracting a regularized control-state reconstruction. Agreement between two disjoint fingerprint gene halves shrinks noisy cell response-strength estimates toward one. Continuous strengths are bounded, normalized to mean one, and regressed on four state factors. Ridge state-by-strength interactions describe coordinated residual gene changes. Source-family estimates are pooled using n/(n+200) reliability; missing/unsupported targets abstain.

At emission, a conditional log2-fold residual is bounded at +/-0.125 after scaling. All population models retain the V16 predicted aggregate mean/pooled targets and its .25 expected + .65 Poisson + .10 empirical mixture before integer moment fitting. Libraries are exact, mean/bulk tolerances are .03/.031 L1. The simple amplitude comparator intentionally changes the aggregate effects by 1.05; it is diagnostic only. Scalar, state-independent and shuffled-state controls share the same learned inputs and count machinery. Zero effects and unsupported targets return V16 exactly.

This is a small statistical factor pilot, not a complete neural generative model or proof of causal state identification. States inferred from treated RNA can retain treatment effects; stable-feature selection reduces but does not eliminate that confounding. Guide efficiency, cell sensitivity and survival cannot be fully separated from these observations. Guide IDs are preserved in original source provenance; this version does not fit a separate guide-efficiency parameter. Continuous scores are response proxies, not calibrated physical knockdown fractions. The query's unobserved responder fraction cannot be recovered from controls alone.

## Tests and reserve

Development: 14 old public target/family pairs, two seeds, 400 genuine treatment cells and predictions, seven fixed arms. The 196 comparisons are repeated draws from three previously studied raw source families, not 196 independent biological experiments. Use unchanged competition-style CPU metric code, family-balanced guards, per-family diagnostics, DE eligibility and paired expression error numerators. Scaled leaderboard scores are not estimated from these raw metrics.

A single conditional finalist must pass the old available-metric guards, improve at least three metrics including DE, reduce centered response-projection Wasserstein distance, and beat scalar, independent and shuffled-state controls on at least two DE metrics each. Exact-null tests must match V16 bitwise. The same rule applies once to the new CD4 raw-cell reserve; missing NMAE is not imputed.

The author-published CIPHER release [Zenodo 21729034](https://doi.org/10.5281/zenodo.21729034) supplies Marson2025_D1_Rest_filtered.h5ad, 1,526,580,299 bytes, MD5 70a49f15fddf62ab767188df0df112de. Verify raw counts and gene IDs before use. Reserve targets are selected alphabetically from metadata with >=400 treatment cells and matched controls, with disjoint within-batch NTC halves. All CD4 evidence, including legacy summary effects, is excluded from reserve-model fitting. This is new raw-cell validation from a study already represented in the source bank, not an independent new study. Reserve treatment scores remain unopened unless development passes.

## Reproduce from Windows CMD

Use the saved frozen plans; do not regenerate them over existing results. The Stack local Python runtime has the required scientific packages; the competition environment is used for official packaging. Metadata string inference is disabled in runner code to avoid the known pandas/Arrow native crash; numerical arrays are unaffected.

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe collect_v22_cd4.py --download
.\.competition_env\Scripts\python.exe test_v22_state_response.py
.\.competition_env\Scripts\python.exe prepare_v22_experiment.py
.\.competition_env\Scripts\python.exe run_v22_experiment.py oracle --workers 4
.\.competition_env\Scripts\python.exe run_v22_experiment.py models
.\.competition_env\Scripts\python.exe run_v22_experiment.py null --workers 4
.\.competition_env\Scripts\python.exe run_v22_experiment.py development --workers 4
.\.competition_env\Scripts\python.exe assess_v22_experiment.py development
.\.competition_env\Scripts\python.exe extract_v22_cd4.py
```

Only after development passes, run confirmation and its once-only assessment. Release scripts must verify the complete model/input/case evidence, count invariants, deterministic replays and official VCC prep. A failed new model retains V16. If an upload-ready fallback is requested, its metadata must clearly identify a byte-identical champion copy; it must not be presented as an improved new predictor.

## Execution notes

Actual scientific runtime is Python3.12.10 / NumPy2.5.4 / SciPy1.18.1 / pandas3.0.6, preserved in V23 analysis/runtime.json and requirements-v23-science.txt. Use the VCC_SCI_PY CMD variable in V23_STRATEGY.md for scientific steps above; the competition environment is reserved for official packaging. V22 failed development and did not open its reserve. Metadata-only source-power probing replaced the underpowered D1 candidate with D2 Rest; the prospective source amendment explicitly uses lane_id, and D1 extraction is archived. See V22_RESULTS.md.
