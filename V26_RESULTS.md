# V26 results: neural count decoder

**No V26 package is released.** The conditional neural decoder fails the frozen H1 real-control calibration gate. Training and count generation completed for the executed pilot, but treatment outcomes and the donor reserve remain unopened. V24 stays the published reference at0.14299785521745675. This is a negative result for this source-trained decoder configuration, not a test of every neural count model.

## What actually ran

Plan SHA256:`142ea7bfeea318ef0f14f94df294a20759b907999c61a7c4f305db9841f96b85`. Fitted 3 prespecified neural models, each for12 epochs, using71,634 sampled source cells (64,922 fitting/6,712 validation) across4 training contexts, with17,472 learned output coordinates. The18533-gene head has explicit masks; unsupported coordinates retain empirical query variation. The complete held-out family is excluded from training, state feature selection, the response basis and the V24 mean predictor. Source data and original cell counts were independently checked.

The neural likelihood models16 latent dimensions and eight state dimensions, with a64/128-unit decoder and16 intervention descriptors. Model epochs are selected by held-cell validation ELBO, which is not a competition metric. The independent, conditional and training-state-shuffled arms share samples and architecture. All fitted weights, source identities, training-only group means and epoch traces are retained.

Executed 6 paired count populations across1 frozen null cases. Each has400 actual controls,400 independent reference controls and400 distinct template controls drawn from the separate inference pool. Six population arms are evaluated. Target-descriptor shuffling is identical to conditional for a zero intervention and reuses that count population; this null screen does not establish target-specific treatment accuracy.

## Measured raw null diagnostics

DE calls use the existing competition-style CPU Wilcoxon/BH scorer. Detection and log-variance errors compare each generated population with the same real reference controls. These are raw diagnostic values, not anchor-scaled leaderboard scores.

|Context|Arm|Generated DE calls|Actual null DE calls|Detection MAE|Log-variance MAE|Mean CPM probability L1 vs V24|
|---|---|---:|---:|---:|---:|---:|
|H1|v24|821|0|0.060798|0.894155|0.00000000|
|H1|empirical|6764|0|0.159986|0.951327|0.00000092|
|H1|independent|389|0|0.068124|0.677430|0.00000095|
|H1|conditional|352|0|0.062203|0.619764|0.00000095|
|H1|target_shuffled|352|0|0.062203|0.619764|0.00000095|
|H1|state_shuffled|337|0|0.067506|0.674665|0.00000092|

Failed checks: **H1_DE_calls, H1_detection, H1_variance**. The mean contract, DE-call envelope, detection envelope and variance envelope are separate checks; see the exact Boolean outcomes in null_H1_result.json. Returning V24 exactly is not treated as calibrated null performance.

The V26 conditional pipeline reduces null DE calls from821 to352 (57.1%) and improves log-variance error from0.8942 to0.6198. Detection error slightly worsens from0.06080 to0.06220; independent real-control detection error is0.01599. Three calibration gates fail; the mean contract passes at about0.000000952 L1 drift. These are one-case raw diagnostics, not treatment gains or a leaderboard prediction. The state-shuffled decoder has337 calls, so the call-count reduction alone does not establish that correct state conditioning is responsible.

The empirical/resampling comparator produces6764 calls and detection error0.15999. Its shared sampling/repair pipeline substantially changes observed-cell distributions. This warns against treating re-sampled observed count probabilities as a neutral empirical reference: the input counts already contain measurement noise. It does not isolate the effects of sampling, stratum constraints and swaps. New generator improvements are collective pipeline effects until separate ablations identify them.

The first numerical sampler designs were rejected on synthetic full-axis mean-contract tests before biological scoring. The executed sampler combines multinomial draws,16 depth-stratum gene-total repairs and bounded integer swaps. All output cell libraries and reference gene totals are exact; mean per-cell expression drift is measured, not assumed zero. Repairs can alter shape, so diagnostics use the final integer counts. A successful mean contract verifies isolation of the mean predictor; it does not establish biological realism.

## Independent audit

- Four focused neural/count tests passed before the plan was frozen.
- 4 source contexts and128 complete original cell rows were checked on the full gene axis.
- 662 fitting-only group means were independently recomputed; validation cells never supplied these means.
- Training/validation identities, state-shuffle batch/partition membership, family exclusions, checkpoint hashes and finite epoch traces passed.
- All executed null-cache hashes, control labels and disjoint actual/reference/inference identities passed.
- 2 full-axis baseline/conditional count replays and their scorer call counts matched exactly.
- Frozen calibration checks were independently recomputed; the automatic stop follows the saved failed gate.

No treatment comparison, later-family training, reserve score or host score is claimed. The pipeline stops at the first failed required calibration as specified before fitting. Existing unexecuted code paths are documented as prepared, not validated. The remaining study/family limits, observational state proxies and source-to-query distribution shift remain relevant.

## Where the remaining errors occur

A separate exact replay of the conditional null locates its352 called genes:

|Coordinate category|All genes|Eligible genes|Null DE calls|
|---|---:|---:|---:|
|Trained; >=50 source detections|15215|10209|142|
|Trained; <50 source detections|2257|89|4|
|Panel genes with empirical fallback|300|292|109|
|Other empirical fallback|761|191|97|

Thus206/352 calls (58.5%) lie in1061/18533 fallback coordinates (5.7% of the output axis). The trained coordinates still contribute146 calls, so changing the fallback alone is not demonstrated to pass calibration. These categories involve different genes; this is descriptive localization, not an isolated causal ablation proving which operation caused each error. The null-gene diagnostic matches the saved conditional count digest exactly. Per-gene q-values, LFCs, detection, variance and source support are in H1_null_called_genes.csv; category counts are in H1_null_error_categories.csv.

This creates two concrete next hypotheses: avoid adding another sampling layer to already observed fallback counts, and adapt learned control distributions using the query context's own available NTCs. Panel genes need exclusion from leakage-prone state features; masking all300 from learned output variation is a conservative choice tested here, not a competition requirement. Any revised mask, sampler or adaptation procedure needs a new plan and validation, rather than changing this failed V26 after seeing outcomes.

## Reporting recovery and latent check

A Path/string mismatch interrupted saving the first completed assessment. The operational assess_v26_recovery.py adapter coerces hash paths, executes the unchanged hash-pinned assessor and preserves table hashes. Original code, fitted arrays, populations and thresholds remain unchanged. run_v26_with_recovery.py is the from-scratch entry point; after recovery the frozen pipeline resumed verified caches and recorded its planned stop. The independent auditor separately recomputes every failed check.

Descriptive source-control diagnostics use1276 validation cells. Conditional aggregated posterior variance ranges0.919..1.062 around the normal generation prior's1, with mean-vector norm0.102. This does not indicate a large diagonal source-prior variance collapse; it does not verify full joint distributions or H1 transfer. See H1_latent_diagnostic.json and root diagnose_v26_latents.py. No diagnostic changes the release rule.

## Implication

This experiment directly tests a complete trainable count-distribution mechanism instead of adding bounded residuals to the old generator. Better source likelihood cannot substitute for calibrated unseen-context distributions. This configuration has not earned a treatment or leaderboard release. Any adaptation to the query controls, change to the likelihood/latent prior or mean-repair rule must be a separately frozen experiment, preserving this result. The present evidence does not identify a hidden competition failure or promise improvement from another architecture. A justified next experiment would calibrate the generator to each query context's available NTC distributions, then learn source-supported treatment-induced distribution changes. This needs a new frozen protocol and matched controls.

## Reproduce

Use the unchanged scientific Python3.12 runtime and versions in analysis/runtime.json. Existing frozen plans refuse overwrite. Source/weight/count arrays are regenerable and ignored by Git. Root scripts and small evidence are versioned. Run:

```cmd
cd /d L:\promec\Animesh\Download\Virtual
set "VCC_SCI_PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
"%VCC_SCI_PY%" audit_v26_experiment.py
"%VCC_SCI_PY%" diagnose_v26_latents.py --family H1
"%VCC_SCI_PY%" diagnose_v26_null_genes.py
"%VCC_SCI_PY%" document_v26_results.py
```

From scratch, V26_STRATEGY.md lists preparation and gated training commands. Replaying the saved pipeline resumes identical fitted weights/caches and preserves its stop. It does not revise parameters or reopen treatment outcomes.
