# V88: converged mixture priors and atom-aware response intervals

The user requested another upload-ready experimental test after V87. This follow-up addresses specific unfinished numerical issues in V86, which improved five of six local D3 endpoints but did not pass its stronger public/shuffled controls. V86/V87 source and evidence are preserved. V66/V82 is the best published reference, overall 0.15088218754984492.

## What changes

The response evidence, covariance family and sampling-error model are inherited from V86. Its 100-iteration EM fits had not converged. V88 fits exactly the same 20,000 source coordinates using a positive-simplex damped Newton solver for the strictly concave penalized mixture log likelihood. The relative KKT tolerance is 1e-8; fitting refuses an unconverged result. Source objectives and sample identities are checked against V86.

The previous normal-moment interval is replaced with numerical quantiles of the actual conditional Gaussian mixture. The null component is a point mass at zero and is retained explicitly in the generalized-inverse CDF. Independent Brent roots check those quantiles. They are plug-in intervals under the fitted model; no calibrated biological coverage is claimed.

The sole production candidate is **exact125**: project the V66 assay response into the mixture's central interval, invert to LFC, then apply one-eighth of the difference clipped at +/-0.5 LFC. The proposal can therefore change each eligible coordinate by at most0.0625 LFC before the existing 80% identity-margin guard. The conservative strength was fixed before new outcomes. Shared125 uses an independently converged shared-only prior; exact125_shuffled fits and conditions on independently shifted source target identities. Neither control can be promoted after outcomes.

Unrolled float64 Cholesky calculations accelerate small Gaussian systems without changing their definition. Fixtures and actual source coordinates compare them against generic NumPy/SciPy solves. Both low-rank context-sensitive and shared covariance components retain an unobserved-context variance term.

## Separation and information limits

Whole recipient studies are excluded before fitting and conditioning: H1; Orion with both lines; CD4 with all donors/states. Query controls supply the allowed baseline embedding and expression scale. Query treatments never enter predictions. Every source/output mask is retained; absent measurements abstain. Proposed outputs exclude panel genes; the existing compositional guard can adjust other final non-own factors. Explicit own-target suppression remains once in the unchanged count emitter.

The public 655-pair screen and D3 donor have already been reused. D3 informed the motivation for this follow-up through V86, so it is not blind validation. Both count seeds use the same donor biology. Large biological sample counts do not turn three independent source study groups into hundreds of independent contexts.

## Prospective decision and upload status

Exact thresholds and hashes are in `experiments/v88_mixture_quantiles/analysis/plan.json`, frozen before fits and outcomes. Strong support requires the original expression/PDS/DE noninferiority checks and superiority to both controls, with both technical seeds. A separate bounded experimental-upload rule was specified prospectively because the user explicitly requested a file to test: small broad-profile/PDS changes, no material donor expression/PDS/FID/JAC/NMAE regression on either seed, and at least two improved endpoints per seed. Failure of strong support remains explicit even if this narrower diagnostic release is permitted. Another arm or strength cannot be promoted after failure.

The count decoder stays V66's residual10. Each generated group has400 cells with original control-template identities and exact original library totals. All360,000 cells cover18,533 genes in competition order. Packaging uses the official VCC prep, with lossless float32 integer-count staging and int64 CSR indices to control memory. Independent validation compares every packaged CSR value/index, verifies every source identity/library, and regenerates12 full400-cell groups.

These are model/integrity tests. Only a new hosted score can establish improvement over V66 on the competition panel. No automatic upload or Git push occurs.

## Root scripts and reproduction

Use the existing scientific Python3.12 runtime at `C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe -X utf8`. Original data/large arrays remain local, with hashes in the protocols. Freeze/fit/output stages refuse to overwrite prior scientific files.

1. `test_v88_mixture.py`; `run_v88_experiment.py freeze`.
2. `run_v88_experiment.py fit --outer H1`, Orion and CD4; repeat each with `--shuffled`.
3. `audit_v88_models.py`.
4. `run_v88_experiment.py prepare`, or equivalently `prepare_v88_context.py NAME` once for each of the five public contexts and CD4_D3_Rest. The latter only partitions the original unchanged calculation for parallel execution.
5. `run_v88_experiment.py means`; `audit_v88_profiles.py predictions` and `metrics`.
6. `run_v88_counts.py counts --offset 0` and620000; optional disjoint `--chunk K --chunks N` partitions. Then `assess --offset OFFSET` and `audit_v88_counts.py --offset OFFSET` for both seeds.
7. `decide_v88_release.py`. Continue only if the prospective experimental-release rule permits it.
8. `fit_v88_submission.py freeze`; `run_v88_experiment.py fit --outer production`; `audit_v88_models.py --production`; `fit_v88_submission.py fit`.
9. `build_v88_submission.py --preflight`, then `build_v88_submission.py --workers 4`.
10. Use the official `.competition_env\Scripts\python.exe -X utf8 package_v88_submission.py`; then the scientific runtime for `validate_v88_submission.py`.

Quantitative results, control comparisons and release limitations are recorded in V88_RESULTS.md when complete. Reproduction requires a separate untouched output folder/checkout; immutable records are not overwritten. Readiness requires both `submission_ready.json` and a passing `independent_validation.json`.
