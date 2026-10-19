# V65 strategy: review the publication history, then test normalized source responses

The reproducible publication inventory is `experiments/v65_closed_response/analysis/published_history.csv`. It contains all 18 saved published versions, including the original root-level receipt and V50's separately named published_result.json. Scores use the same panel and anchors. Ranks are historical snapshots affected by other entrants; compare scores to assess our own improvements. Changes between versions often combine components and do not prove individual causal effects.

V10's pooled response/population change and V11's consensus coverage produce the largest later gains. V13/V16/V24 refine that behavior; V63 provides the first material improvement after the V24 plateau, primarily in expression accuracy. Stacking (V20), retrieval correction (V36), response priors (V50), fixed blending (V55) and state transfer (V61) fail to improve the best prior published score. V24 improves all six components over V16; V63 improves four over V24 while losing fidelity and reach.

## Hypothesis fixed before scoring

V63 uses source profiles covering about 68%-100% of a CPM library without reclosing the available gene axis. It converts LFC to assay-scale deltas without explicit source perturbation closure. V65 tests whether those modeling choices leave avoidable cross-assay offsets. This is a hypothesis about a changed estimator, not proof of a bug: incomplete axes do not recover missing biology by normalization.

Three candidates isolate available-axis reclosure, reclosure plus source perturbation closure, and a stronger correction of the latter. Each uses the inherited source uncertainty/shrinkage, source-only leave-one-target-out common-response removal, query inversion, and the original V24 margin guard. The own-target factor is 0.1 inside source closure and exactly once in the query emitter; stored response effects keep own coordinates zero. Diagonal uncertainty ignores cross-gene normalization covariance. The unchanged V24 residual10 generator allows the response-estimator comparison without a simultaneous dispersion change.

Select the smallest equal-study target-weighted desired-profile error ratio to V63 subject to the frozen mean/PDS guards. Lock before any new count results; no post-count alternate promotion. Compare the selected candidate, V63, V24 and a cyclic source-target derangement using common templates and seeds. Exact old baseline counts may be reused only with verified receipt/array hashes. Use the first case for each unique target/context: 30 pairs, five contexts and three studies, two technical generation seeds. Reused public biology is not independent confirmation. D3 remains unopened.

The official metric uses log1p profiles of group-summed counts normalized to 50,000, excludes the intervention's own gene from scoring, and removes all 300 panel targets for PDS. [Official metric specification](https://raw.githubusercontent.com/ArcInstitute/cell-eval2/main/docs/vcc2026_metrics/vcc2026-metrics.md). The current website evaluation route provides only its client-side loading shell to the text reader. Local raw diagnostics have no access to the hidden official anchor scales.

## Reproduction

Use the existing scientific environment and requirements-v35-science.txt; `.competition_env` is for official packaging. `SCIENCE_PYTHON` denotes the documented Python 3.12 scientific interpreter. Freeze once in a fresh directory; do not overwrite prior plans. Source/cache prerequisites are the existing V7/V24/V34/V37/V63 documented workflows.

```cmd
SCIENCE_PYTHON review_v65_history.py
SCIENCE_PYTHON test_v65_response.py
SCIENCE_PYTHON v65_closed_response.py freeze
SCIENCE_PYTHON v65_closed_response.py prepare
SCIENCE_PYTHON audit_v65_boundaries.py
SCIENCE_PYTHON run_v65_experiment.py means
REM Only if the frozen response screen selects a candidate:
SCIENCE_PYTHON fit_v65_submission.py freeze
SCIENCE_PYTHON run_v65_experiment.py parallel
REM After complete count assessment, fit normally if supported;
REM --experimental is required otherwise and records failed guards.
SCIENCE_PYTHON fit_v65_submission.py fit --experimental
SCIENCE_PYTHON build_v65_submission.py --preflight
SCIENCE_PYTHON build_v65_submission.py --workers 4
.\.competition_env\Scripts\python.exe package_v65_submission.py
SCIENCE_PYTHON validate_v65_submission.py
```

`prepare_v65_release_tools.py` generated the root writer/package/payload-validator adapters from reviewed V63 files without changing predecessors. `verify_v65_effects.py` replays every production response and checks all prior-margin inequalities; the payload validator checks all packaged values/axes, original control provenance/libraries and 12 full-group replays. Scientific files and compact receipts are locally committed; large data, arrays and packages are ignored. No automatic upload or Git push.
