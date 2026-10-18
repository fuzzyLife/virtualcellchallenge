# V60: control-only adjustment of the state representation

V60 preserves V59 as an independent experiment. Its only scientific change is the cell-state normalization: fit log1p(CPM) on log nonpanel raw depth, log nonpanel detected-gene count, their squares and interaction, using NTC cells within each lane. Use ridge .1 with unpenalized intercept; standardize residuals by NTC SD with floor .5 and clip +/-5. Both nuisance covariates exclude all 300 panel target genes. All fitting applies to source controls only; treated measurements never estimate these coefficients.

The motivation is the NTC-only V59 state-quality audit: CD4 PC1/PC2 correlate .87-.92 with library depth and .89-.94 with detected genes. This does not prove those axes are technical. Biological cell size and transcriptome complexity may be removed too; V60 must be compared experimentally rather than called a guaranteed correction.

The 512 stable genes, eight PCA coordinates, four states, whole-study exclusions, two raw views, density/conditional estimators, same-target donor weights, paired shuffle, 50% correction and release gates are unchanged. V60 has separate state models, mechanisms, predictions and result directories. Verified V59 raw CPM caches and V52 row identities are reused without re-extraction. The CPM closure and missing-gene routing fixes are inherited. All seeds remain paired with V59. No diagnostic arm can be promoted after outcomes.

V60 primary joint50 must pass both its response and integer-count tests; it is not selected because the new representation looks attractive. D3 remains unopened. No automatic upload or Git push.

## Reproduce

```bat
cd /d L:\promec\Animesh\Download\Virtual
set "SCIENCE_PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
"%SCIENCE_PY%" -X utf8 test_v60_states.py
"%SCIENCE_PY%" -X utf8 run_v60_experiment.py response
"%SCIENCE_PY%" -X utf8 run_v60_experiment.py counts
"%SCIENCE_PY%" -X utf8 v60_depth_adjusted_states.py plot
```

First generate V59 raw caches using its corrected controller. Plans refuse overwrite. Numerical contracts demonstrate that changing treatment data cannot alter the fitted NTC regression and that a synthetic independent signal survives depth adjustment. Those unit checks are not evidence of biological predictive improvement.

## Observed stopping scope

All 366 response cases completed, but the matched shuffled/pooled control criterion failed in every study. No V60 count experiment or improved submission is claimed. Population-level reconstruction diagnostics are separately recorded as a post-hoc correction of the original tiny-lane error aggregation; no V60 predictor or threshold was changed. V61 isolates pooled state estimation using the exact V60 state basis.
