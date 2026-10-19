# V69: probability-preserving count bridge

Frozen before source fitting and count outcomes. The independent implementation uses a Gamma-Poisson posterior coupling, with beta/gamma shape changes, binomial removals and Poisson additions. It tests the distribution-learning hypothesis in [the research review](POST_V68_RESEARCH.md), not an implementation of PopPert or ExpressionVAE. Zero response exactly preserves original counts. No exact row-total or gene-total repair follows. The working NB marginal assumption and depth-conditioning are approximations; empirical control counts need not follow this model.

The mean predictor remains frozen V66. Five public raw contexts supply target/gene dispersion shifts, split-view sign support and study-balanced shrinkage. Whole-study exclusions remove both Orion lines or all CD4 conditions. Panel genes do not provide dispersion labels. Production uses only public source treatments and A/B/C controls. Own suppression is applied once.

A labelled-data mean/decoder diagnostic precedes deployable testing. It uses ten cases, splitting treatment populations into separate fitting/scoring halves (100 or200 cells), five diagnostic arms, and explicitly reports realized means. These oracles are not production predictors. The deployable screen compares mean-only bridge, half-strength learned shape bridge and a population with25% bridge cells/75% unchanged V66 cells, plus shuffled-shape evidence. First two unique targets/context select the arm; remaining unique targets confirm it. Reused biology is not a blind holdout. The hybrid cannot inherit the full bridge's exact null identity.

All formulas, seeds, cases, selection rules and strict confirmation guards are stored in `experiments/v69_count_bridge/analysis/plan.json`. Six numerical tests check count identity, reproducibility, integer values, marginal simulation, unmeasured-zero abstention and own suppression. Source moments are method-of-moments estimates, not a fitted NB likelihood or identified causal intervention.

The user requested an upload artifact. If no arm meets the development eligibility checks, the predeclared hybrid25 fallback is labelled experimental; confirmation does not trigger post-hoc selection or a claim of improvement. All failed guards and oracle limitations must be reported with the release. No automatic upload occurs. The untouched donor reserve stays unopened.

Commands use the pinned scientific Python (`%SCIENCE_PY%`):

```cmd
"%SCIENCE_PY%" -X utf8 test_v69_count_bridge.py
"%SCIENCE_PY%" -X utf8 run_v69_experiment.py freeze
"%SCIENCE_PY%" -X utf8 run_v69_experiment.py fit
"%SCIENCE_PY%" -X utf8 run_v69_experiment.py diagnostic
"%SCIENCE_PY%" -X utf8 run_v69_experiment.py assess-diagnostic
"%SCIENCE_PY%" -X utf8 run_v69_experiment.py development
"%SCIENCE_PY%" -X utf8 run_v69_experiment.py select
"%SCIENCE_PY%" -X utf8 run_v69_experiment.py confirmation
"%SCIENCE_PY%" -X utf8 run_v69_experiment.py assess-confirmation
"%SCIENCE_PY%" -X utf8 run_v69_experiment.py null
```

Source counts, diagnostics and model arrays are regenerable ignored files; compact receipts, scripts and documentation are committed. No new leaderboard result is claimed.

## Executed additions and production runtime

The original scientific runner above remains byte-frozen. After the labelled diagnostic exposed undefined NMAE, `assess_v69_confirmation.py` was added to preserve unavailable values and fail their guards explicitly. Before development, `run_v69_matched_control.py` was frozen to add a quarter-mixture shuffled control. The original full shuffled population is not a matched comparator for hybrid25. The separate representative-control null plan preserves and explains the prefix-template sampling results.

Actual deployable evaluation uses content-addressed memoization (`run_v69_cached.py`) through the bounded two-process `run_v69_parallel.py`, then the matched control and safe assessor. Memoization keys contain complete CSR values, indices, pointers, shape/dtype and gene axis; exact-statistic and mutation-key tests pass. The confirmation control may wait for completed input receipts via `run_v69_matched_wait.py`; no missing values are fabricated.

The experimental fallback was locked before confirmation. Production parameters were fitted/generated while confirmation ran, with assessment explicitly pending and official packaging blocked until finalization. `finalize_v69_confirmation.py` updates assessment metadata only after all count generation and confirmation finish, preserving prior records and proving unchanged model-array bytes.

Production hybrid25 samples the selected100 bridge rows directly. All Gamma/Beta/Binomial/Poisson transitions act independently across rows conditional on fixed per-row original library and external p/alpha/effect/dispersion parameters. Subsetting preserves their joint probability law with a new documented random-stream realization. The earlier400-row implementation, metadata and unused partial output are retained; source parameters, candidate and validation thresholds are unchanged. New preflight and deterministic regeneration check the actual deployed sampler. The frozen public evaluations retain their original400-row draws.

```cmd
"%SCIENCE_PY%" -X utf8 run_v69_parallel.py development
"%SCIENCE_PY%" -X utf8 run_v69_matched_control.py development
"%SCIENCE_PY%" -X utf8 run_v69_experiment.py select
"%SCIENCE_PY%" -X utf8 run_v69_parallel.py confirmation
"%SCIENCE_PY%" -X utf8 run_v69_matched_wait.py confirmation
"%SCIENCE_PY%" -X utf8 assess_v69_confirmation.py
"%SCIENCE_PY%" -X utf8 build_v69_submission.py fit
"%SCIENCE_PY%" -X utf8 verify_v69_model.py
"%SCIENCE_PY%" -X utf8 build_v69_submission.py preflight
"%SCIENCE_PY%" -X utf8 build_v69_submission.py build
"%SCIENCE_PY%" -X utf8 finalize_v69_confirmation.py
"%SCIENCE_PY%" -X utf8 audit_v69_release.py results
.\.competition_env\Scripts\python.exe package_v69_submission.py
"%SCIENCE_PY%" -X utf8 audit_v69_release.py payload
"%SCIENCE_PY%" -X utf8 document_v69_release.py
"%SCIENCE_PY%" -X utf8 seal_v69_release.py
"%SCIENCE_PY%" -X utf8 seal_v69_release.py --check
```

These stage commands refuse overwriting completed plans/models/artifacts. For the recorded parallel execution, production `fit` used `--allow-pending-confirmation` only after hybrid25 was locked as fallback. A sequential reproduction after confirmation does not need that flag. Historical code and addendum receipts live in `analysis/execution_history` and the corresponding named JSON records; the unused incomplete H5AD is ignored under `runtime_history`.
