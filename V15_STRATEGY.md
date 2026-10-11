# V15: nested response ensemble and source-linked population transport

V13 remains the published reference. V15 tests response prediction and population generation separately, combines only predeclared candidates, and performs reserved algorithm verification before packaging. No upload is performed by the scripts.

## Response experts and ensemble

`v15_ensemble_methods.py` implements three experts: the exact V13 unknown-knockdown pooled response, a Student-t reliability variant with the same unknown-knockdown handling, and a similarity-weighted source specialist. Source reliability uses transfer errors between training families; the specialist uses pairwise-valid Spearman correlations of unperturbed profiles, excluding all 300 panel genes. It does not assign literal cell-line identities from marker correlations.

A six-point convex-weight grid retains at least 50% V13. For each outer-held family, weights are selected only from inner-family held-out predictions within the remaining sources. All transformations, source-reliability estimates and response fitting exclude the appropriate outer and inner held families. CD4 conditions average within one source family; they are not extra independent votes. Reserved target names are omitted from inner summary selection and the outer summary comparison, although their measurements in other training contexts remain legitimate known-target inputs. Training-source reliability uses all available known-target measurements; this is not unseen-target validation.

Source measurements and invalid entries are masked explicitly; missing values abstain. Effects are bounded to +/-3 log2 and own-target downstream effects remain zero. Direct own-target repression is handled by the emitter as in V13. Error correlations and cases where each expert corrects/breaks V13's reference-DE directions are reported. These diagnostics do not establish biological correctness on private A/B/C responses.

The screen freezes either the nested ensemble or the fixed 50% Student/V13 blend. Summary objectives are proxies; actual generated counts determine adoption. The Student expert here is not a bitwise copy of the older Student3 candidate: unquantifiable source knockdown now receives V13's neutral handling.

## Population model

`v15_population_methods.py` extracts five relative-expression quantiles in two control-defined depth groups from public H1, HCT116 and HEK293T treatments. Each source population is divided by its own mean proportions, keeping predicted treatment averages separate from relative distribution changes. Controls use at most 800 distinct cells per depth group. Treatments need >=80 unique cells and >=25 in each depth group; genes must be measured and exceed 5 CPM in both source control and treatment. Insufficient support abstains.

Quantile differences are capped at +/-4 relative units. Distribution pooling excludes the scored family, uses bounded source response confidence and cell-count reliability, and removes own-target shape changes. K562/CD4 currently supply response summaries but no marginal-distribution arrays in this experiment.

For a query, sampled empirical controls are divided by their context-wide mean proportions. Within low/high template-depth groups, ties receive reproducible randomized ranks; interpolated source quantile differences are added with 0.5 strength and clipped at zero. This is a bounded rank-dependent marginal adjustment, not an optimal-transport solver or a guarantee of preserved copula/order after the adjustment. The result is scaled to predicted treatment means and blended into the compact population only on measured source-supported genes, at 10% or 25%. Dual-moment integer fitting then matches mean-CPM/pooled proportions while preserving each sampled library exactly. Final fitting can change the initially proposed marginal distribution, which is why generated-count tests matter.

This is distinct from uniformly restoring control variance, which failed earlier. It still risks transferring sampling differences and source-batch/assay effects. Extreme quantiles are noisy and sensitive to unequal source population sizes; bounded changes, modest blends, negative controls and held-family evaluation are safeguards rather than proof of transferability. Source quantiles are pooled across source batches and stratified by depth; this does not remove residual batch confounding. No held-family treatment data or hidden competition treatments enter deployable fitting.

## Frozen validation

`experiments/v15_weighted_transport/analysis/plan.json` fixes:

- Development: 14 targets (H1 4, HCT116 4, HEK293T 6), two seeds, 400 genuine treatments and 400 predictions per target/seed.
- Six count arms: exact V13, selected response/compact, V13/quantile25, selected response/quantile10, selected response/quantile25, and selected response with shuffled-target quantile25.
- Development therefore generates 168 populations. No duplicated treatments are used to claim added sample size.
- Nine reserved verification target names across all three families; one selected candidate versus V13 at two seeds gives 36 populations if a candidate advances.
- Null: 18 populations from disjoint 200-cell template/actual NTC draws, with independent right controls. Null query means use the full left pool, including those NTC actual cells; the selection check uses generated-versus-independent-right DE calls, not the optimistic null-truth expression diagnostic.
- All primary references use full independent right controls. V14 control halves, actual treatment identities and seeds are retained, but the saved control-pool order changes template sampling; the original baseline is independently replayed on V15 inputs (see the audit correction below); old public biology is reused and is not blind confirmation.

The response screen is followed by generated-count development. Require improvement in at least three available metrics, including fidelity/reach/overlap, while guarding discrimination, fold changes, expression and worst-family DE declines. Missing NMAE panels are explicitly reported; they prevent a claim of complete six-metric verification, but do not fabricate zero errors. Tiny expression denominators are protected using paired absolute error numerators. The prospective 2% tolerance is exploratory, not statistical noninferiority or a retroactive V14 rescue. Passing candidates are ranked by equal metric-improvement ranks; no private anchored score is estimated.

The original shuffled25 arm is exactly matched to response/quantile25 but differs in blend from response/quantile10 and in response means from V13/quantile25. Before any generated-count outcomes, `matched_control_plan.json` additionally requires an exact response-mean/blend-matched shuffled-target check for any distribution finalist, winning at least two of fidelity/reach/overlap. This extra check can reject a finalist but cannot retune or reselect one. It preserves the original frozen screen/count code and records its timing and hashes separately.

A once-only assessment freezes the candidate before reserved verification. Public source biology has already been inspected in earlier versions; disjoint current target names and whole-family model fitting do not make these new biological experiments. No hidden leaderboard improvement is promised. A candidate with incomplete primary metrics can only be described as exploratory. Failed verification retains V13.

## Reproduction from Command Prompt

Preserve existing frozen plans and outputs; preparation and once-only assessment refuse to overwrite them. Source/model/count caches resume only when their input and byte hashes match.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m unittest test_v15_ensemble -v
.\.competition_env\Scripts\python.exe prepare_v15_experiment.py
.\.competition_env\Scripts\python.exe run_v15_experiment.py sources
.\.competition_env\Scripts\python.exe run_v15_experiment.py screen
.\.competition_env\Scripts\python.exe run_v15_experiment.py null --workers 4
.\.competition_env\Scripts\python.exe run_v15_experiment.py development --workers 4
.\.competition_env\Scripts\python.exe assess_v15_experiment.py development
```

Only after a new candidate advances:

```cmd
.\.competition_env\Scripts\python.exe check_v15_matched_control.py
.\.competition_env\Scripts\python.exe run_v15_experiment.py confirmation --workers 4
.\.competition_env\Scripts\python.exe assess_v15_experiment.py confirmation
.\.competition_env\Scripts\python.exe replay_v15_baseline.py
.\.competition_env\Scripts\python.exe audit_v15_experiment.py
.\.competition_env\Scripts\python.exe fit_v15_submission.py
.\.competition_env\Scripts\python.exe build_v15_submission.py --preflight --workers 4
.\.competition_env\Scripts\python.exe build_v15_submission.py --workers 4
.\.competition_env\Scripts\python.exe validate_v15_submission.py
```

The supplemental matched-control freeze is versioned and pins its script; reproducing from saved plans uses that file. An independent experiment audit checks original identities/splits, nested grid selection, all cached profile hashes and manually sorted/interpolated source quantile samples. The release validator checks every packaged CSR value/index/pointer/axis, all original control identities/guides/library totals and twelve deterministic complete-group regenerations. Official VCC prep validates the upload format. These checks establish integrity, not biological accuracy.

Large arrays, counts, packages and complete sampling provenance are excluded from Git. Root scripts and small JSON/CSV evidence are tracked. Numerical V13/V14 code and prior packages remain unchanged. Measured outcomes and the eventual release decision will be recorded in `V15_RESULTS.md` and, if a new package passes, `V15_SUBMISSION.md`.

## Audit correction and observed release outcome

All four deployable development candidates failed the frozen policy. The 50/50 response blend improved eligible NMAE by 0.002245 but worsened fidelity by 0.001270, reach by 0.003969, overlap by 0.000032 and paired expression error by 2.23%. Distribution blends failed the control-only false-positive guard as well. Reserved verification and production fitting were not opened. V13 remains the published/upload reference; resubmitting its identical file is unnecessary. Full measured evidence is in `V15_RESULTS.md`.

The initial audit incorrectly expected 18 historical V14 count hashes to match. V14 sampled unsorted control pools before saving them sorted; V15 samples the saved sorted pools. Actual treatment identities, complete control halves and seeds match, but the sampled templates differ. Within V15, every candidate shares the exact same template, actual and reference cells. The historical mismatch is preserved in `audit_correction.json` and `historical_baseline_differences.csv`. `audit_v15_experiment_initial.py` preserves the failed auditor. `replay_v15_baseline.py` verifies the unchanged original loader, gating and emitter on V15's exact inputs, including count hashes, profiles and DE results; run it before the corrected independent audit even when development rejects all candidates. Frozen numerical methods, the plan, outcomes and the rejection remain unchanged.

After the audit, generate the measured reports and README section with:

```cmd
.\.competition_env\Scripts\python.exe document_v15_submission.py
```

A mechanistic limitation of the tested population blend is that it also multiplies the existing Poisson coefficient by `(1 - alpha)`, reducing 0.65 to 0.585 or 0.4875 on supported genes before moment fitting. This is not an isolated addition of dispersion. The zero-delta NTC arms already increase false calls, so the failure cannot be attributed solely to source-specific quantile differences. A noise-preserving alternative requires a separately frozen experiment; none was silently substituted after the outcomes.
