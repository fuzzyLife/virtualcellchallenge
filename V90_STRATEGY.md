# V90: explicitly post hoc leaderboard diagnostic

The user requested a new submission to test. V90 packages the unchanged V86 `interval25` response method after completing its second count-generation-seed test as V89. It is **not a new architecture or a validated improvement over V66**. It is an explicitly post hoc diagnostic choice from known public results.

V88 and V89 remain rejected under their original frozen research protocols. Their plans, results, failed checks and no-package decisions are preserved. V90 does not change those thresholds or label their failures as successes. A separate diagnostic release was chosen because the candidate improved expression-profile error, direction fidelity and fold-change error on both donor seeds, with unchanged discrimination, while significant-gene overlap was unstable. The user-requested hosted test can measure that aggregate tradeoff on the competition panel; local data cannot establish it.

## Known evidence and limitations

The 655 reused public pairs show0.3734% lower desired-profile error and a raw PDS change of-0.00053622. The candidate did not beat the shared/shuffled controls on every criterion. These studies have informed many project decisions.

The 14-target donor seed0 result is reused from V86; its first-seed outcome informed selection. V89 emitted42 new primary/control populations on seed620000 and reused14 verified V84 baselines. Both seeds use the same donor biology. There is no fresh biological validation, and no new seed0 emission is claimed.

Both seeds improve raw profile MSE, fidelity and NMAE; PDS is unchanged. Jaccard improves in seed0 and decreases in seed620000. The latter caused V89's prospective release gate to fail. V90 deliberately exposes this tradeoff as a diagnostic; it does not assert that a weighted average of these raw numbers predicts the official anchored Overall score.

V88's numerically converged mixture fit failed overlap/discrimination checks. V90 therefore retains the original finite-iteration V86 mixture fit and its normal-moment interval projection. Those approximations are not claimed converged or biologically calibrated. The source data, covariance definitions, quarter-strength correction and80% identity guard are unchanged. Production fits all five measured source families in an isolated Python module namespace; original V86/V89 files are not modified. Float64 Cholesky arithmetic is used only after parity checks on real source coordinates.

## Production and integrity

`v90_diagnostic_release.py freeze` records the known evidence hashes and explicitly labels selection post hoc. It requires the stable three-metric gains and unchanged PDS actually to exist in both saved donor seeds. `fit_v90_submission.py freeze` then hashes the production adapters before fitting.

The production model uses the immutable V66 A/B/C effects as its baseline, original query control profiles, original control-cell selections, original libraries, and the same seed20261011. `v16_methods.emit` residual10 is unchanged, including explicit own-target suppression once. The output has360,000cells,900context-targetgroups,400cells per group, and18,533genes in competition order.

The V86 implementation computes control outputs alongside the primary, so both all-source real/shuffled priors are fitted, but only interval25 is deployed. Public-control study exclusions and numerical model audits come from the preserved source experiments. `audit_v90_production.py` checks all-source Gaussian arithmetic, component weights and actual finite-iteration training histories without claiming convergence.

The package is prepared using the official CLI. Lossless float32 integer-count staging with int64 CSR indices controls memory. Independent validation decompresses the official payload, compares every count/index/pointer with the generated source, verifies every original control identity/library and all axes, and regenerates12 complete400-cell groups exactly. File integrity does not establish biological superiority.

## Reproduce

Use the existing scientific Python3.12 runtime at `C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe -X utf8`, except official packaging uses `.competition_env\Scripts\python.exe`. Original data and large scientific arrays remain local and hash-recorded.

```text
v90_diagnostic_release.py freeze
fit_v90_submission.py freeze
v90_diagnostic_release.py fit-production
v90_diagnostic_release.py fit-production --shuffled
audit_v90_production.py
fit_v90_submission.py fit
build_v90_submission.py --preflight
build_v90_submission.py --workers 4
package_v90_submission.py
validate_v90_submission.py
```

`finish_v90_release.py` executes this sequence and never uploads or pushes. Use an untouched separate output checkout/folder for reproduction; existing plans, models and predictions cannot be overwritten. `report_v90_release.py` writes the final result/README evidence after validation. Source V88/V89 numerical/count experiments remain reproducible through their root strategy documents.

The upload file is `experiments/v90_interval_diagnostic/prediction.vcc` only once `submission_ready.json` and a passing `independent_validation.json` both exist. V66/V82 remains the best published score until a new hosted result is available. All noninferiority failures and uncertainty remain explicit.
