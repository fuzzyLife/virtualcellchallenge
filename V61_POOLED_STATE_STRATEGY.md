# V61: pool matched state evidence before shrinkage

This is a separately frozen estimator following V59/V60. Their original files, scores and rejected candidates remain available. V61 reuses the exact V60 nuisance-adjusted state bases for the public folds. All 300 panel targets remain excluded from state features and nuisance covariates; complete recipient studies remain excluded from basis/response fitting.

The important statistical change is the unit of estimation. V52 metadata allow as few as two treated cells in an individual target/lane view. V59/V60 estimated nonlinear state ratios and shrank conditional effects separately in each lane, before pooling. This nearly erased weak lane-state effects. Their original oracle table also averaged tiny-lane errors; it is not a reliable population reconstruction diagnostic and must not be compared numerically to V61's aggregate errors.

For each target and disjoint view, V61 sums treated counts and expression by state across lanes. Conditional controls remain matched to each treated lane/state, weighted by treated state counts. Abundance controls average NTC state fractions using the target's treated lane counts. Then estimate density ratios and conditional fold changes, and apply the same priors/bounds as V59/V60. Two views combine with treated-state count weights. An empty individual lane/control state uses that lane's overall control mean before pooled support checks; this is an explicit statistical imputation assumption, not an observed counterfactual. Entirely unsupported pooled states and genes abstain.

The corrected oracle fits one full matched population view and predicts the other, computing errors after population aggregation. It cannot enter deployable predictions. Direct mean, zero-response, within-only, composition-only, joint and state-permuted diagnostics are retained. This revised diagnostic was specified after noticing V59/V60's tiny-lane errors; it is not retroactively called the original prespecified oracle.

The transfer formula, same-target evidence, source masks, rank routing, study weights, 50% correction, target shuffle, fixed count cohorts, count emitter and release gates are unchanged. Fixed primary candidate: **V61 joint50**. Diagnostic arms cannot be selected post hoc. D3 remains unopened.

## Release and reproduction

The user requested a new release. The production candidate is fixed before its count outcomes. Public guard results determine whether it is supported as a replacement. `--experimental` explicitly permits packaging the fixed candidate for diagnostic upload after failed gates; such a file must report the failures, preserve V24 and make no improved-score claim. Nothing is automatically uploaded.

```bat
cd /d L:\promec\Animesh\Download\Virtual
set "SCIENCE_PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
"%SCIENCE_PY%" -X utf8 test_v61_pooling.py
"%SCIENCE_PY%" -X utf8 run_v61_experiment.py response
"%SCIENCE_PY%" -X utf8 run_v61_experiment.py counts
"%SCIENCE_PY%" -X utf8 run_v61_release.py --experimental --workers 4
```

The release runner verifies public boundaries, compares an independently implemented production kernel with all saved public predictions, fits on permitted public sources, reconstructs all 900 production effect vectors through the original transfer implementation, checks preflight counts, generates integer cells, runs official VCC prep and compares every packaged sparse value/index/pointer and all cell/gene axes. It verifies all original control identities/libraries and replays 12 complete groups. Unchanged groups reuse the independently validated V24 counts exactly. Large arrays/packages remain outside Git; root code, receipts and reports are committed locally.

## Cache durability and prerequisites

The scientific commands use the project's verified V52 labels, V58 baseline exports, original public raw source matrices and provided competition controls. They are prerequisites, not files downloaded by the V61 runner. Earlier collectors/protocols document their acquisition. Git contains code, compact evidence and reports; large original datasets, model arrays and upload packages are ignored.

`export_v61_arrays.py` copies V61 bases/mechanisms/predictions/count statistics and the reused V59 raw CPM caches into `analysis/arrays/`, with a checksum manifest. `restore_v61_arrays.py` verifies and restores those exports after temporary-directory cleanup without overwriting a different file. Production model arrays and query inputs already live in the project. Reproduction on another machine needs the source-data prerequisites or exported arrays, plus Python3.12.10 and the pinned requirements-v23-science runtime (with its existing plotting dependencies); official prep uses the separate competition environment.

Public population-count tests use moments from the sampled template, whereas the unchanged production generator uses all provided NTCs for its context moments, following V24. This is a material validation/production difference. The limited 4/8-target common-control panel diagnostic is not a substitute for official 300-target PDS or private replicate anchors. All public scores remain descriptive.
