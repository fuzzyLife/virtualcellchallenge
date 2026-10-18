# V62: source-common response removal and assay-space transfer

V61 lost 18.35% versus the V24 overall score. This experiment changes response estimation while retaining the established V24 integer-count generator. It does not reuse V61 state corrections.

## Research and hypothesis

The [Molina and Zhang response-decomposition preprint](https://www.biorxiv.org/content/10.64898/2026.07.24.740459v1) and [author ANOVA implementation](https://github.com/xinyizhanglab/perturbation-decomposition/blob/main/decomposition/anova.py) separate global, context, conserved-target and interaction components. Their findings motivate testing which components can transfer; they do not establish that our source-common component is technical noise or predict a VCC improvement. Our estimator is an independent masked leave-one-target-out adaptation, not their full method. The author's four-context ANOVA describes measured data and cannot supply an unobserved query template.

The [official metric specification](https://raw.githubusercontent.com/ArcInstitute/cell-eval2/main/docs/vcc2026_metrics/vcc2026-metrics.md) forms pooled-count profiles as log1p(50000 * gene proportion). Moving a proportional log fold change between very different expression baselines changes its importance in that space. We therefore compare transferring effects directly in those expression units. Source profiles here are mean per-cell CPM proxies, while final evaluation uses actual pooled counts; this mismatch is explicit.

Three fixed candidates: centered50 removes source-common LFC programs before pooling; assay50 pools source log1p-profile differences and inverts them against query controls; assay_center50 combines both ideas. All use half-strength corrections bounded to 0.5 log2 units from the fold-specific V24 mechanism. Missing source coordinates abstain. Measured zeros remain zeros. The separate established emitter applies own-target knockdown once.

## Validation

The exact plan is experiments/v62_conserved_response/analysis/plan.json. Whole H1, Orion (HCT116 and HEK293T together), and CD4 studies are excluded from each response bank, its support decisions and its templates. Templates exclude the intervention currently being predicted and require 20 other supported targets. No recipient treatment data enter inference.

The same 44 previously examined raw cases and audited V37 populations are used. This is reused public development, not independent new biological validation. Mean-screen selection is locked before count development. Counts compare V24, the selected candidate and cyclically shuffled source-target evidence, with common controls, templates and seeds. Confirmation uses the remaining cases and a second seed only after development support. Small-panel PDS is a diagnostic and is not a forecast of 300-target PDS. Raw local errors are not the CLI's scaled accuracies. Missing NMAE stays missing; no private anchors are estimated. D3 remains unopened.

A pre-score fixture correction replaced exact floating zero with 1e-12 tolerance after expm1/log1p round-trip errors of at most 3.3e-16. Numerical predictions were unchanged; original plan and fixture are archived alongside the final plan.

## Reproduce

Use the existing Python 3.12 scientific environment and the preceding V7/V24/V34/V37 input workflows. Root scripts verify saved input hashes. From CMD, replace SCIENCE_PYTHON with that environment's interpreter:

```cmd
SCIENCE_PYTHON test_v62_response.py
SCIENCE_PYTHON v62_conserved_response.py freeze
SCIENCE_PYTHON v62_conserved_response.py prepare
SCIENCE_PYTHON v62_conserved_response.py means
SCIENCE_PYTHON run_v62_counts.py development
SCIENCE_PYTHON run_v62_counts.py assess-development
SCIENCE_PYTHON run_v62_counts.py confirmation
SCIENCE_PYTHON run_v62_counts.py assess-confirmation
```

Freeze only in a new experiment folder; preserve the original plan. The numerical-tolerance amendment above is already included in this run. No automatic upload or Git push.
