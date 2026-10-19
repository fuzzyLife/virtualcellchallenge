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
