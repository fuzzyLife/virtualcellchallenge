# V27 results: control-anchored relative responses

**No V27 submission is released.** The fixed finalist failed prospective checks. V24 remains the best recorded published submission, overall 0.14299785521745675.

Frozen plan SHA256: `a1e97800b9d6ca3778791c2612b457ac8ee3d29400d057ea17ad8d303c11212f`. Historical V24--V26 methods, packages and numeric evidence remain unchanged.

## Sampling-stage ablation

|stage|null_DE_calls|detection_MAE|log_variance_MAE|mean_reference_L1|
|---|---|---|---|---|
|raw_template|0|0.016739|0.128291|0.012578|
|raw_second_multinomial|3225|0.091642|0.521762|0.020524|
|tilted_multinomial|3176|0.092043|0.522450|0.016430|
|depth_repaired|6750|0.159737|0.946338|0.001544|
|CPM_repaired|6764|0.159986|0.951327|0.000001|
|v24_reference|821|0.060798|0.894155|0.000000|
|raw_anchor_projected|240|0.084583|0.741372|0.000001|

This spent H1 control case uses the same real reference throughout. Stages are sequential ablations, not independent studies. Intermediate stages need not meet V24 mean constraints; their drift is explicit. The final repaired empirical population exactly reproduces the frozen V26 count digest. The raw control matrix already contains measurement noise. Added sampling and count repairs can each alter its distribution; neither is a neutral empirical reference.

The descriptive artifact check also records signed detection errors and support changes:

|stage|signed_detection_error|changed_molecule_fraction_vs_raw|panel_DE_calls|other_DE_calls|
|---|---|---|---|---|
|raw_template|-0.006054|0.000000|0|0|
|raw_second_multinomial|-0.091908|0.131077|75|3150|
|tilted_multinomial|-0.092352|0.131426|82|3094|
|depth_repaired|-0.160058|0.165363|169|6581|
|CPM_repaired|-0.160303|0.170309|171|6593|
|v24_reference|0.060015|0.185049|18|803|
|raw_anchor_projected|-0.084949|0.050443|7|233|

The second multinomial draw removes409,390 nonzero cell-gene entries in this case. Depth-stratum repair subsequently loses324,340 and adds27,220 entries: a net loss297,120. These are measured sequential support changes; the later V24-reference and raw-anchor rows are separate comparators, not causal continuation stages.

## H1 control calibration

|arm|mean_calls|max_calls|detection_MAE|log_variance_MAE|max_mean_L1|
|---|---|---|---|---|---|
|v24|982.333333|1128|0.059941|0.872152|0.000000|
|v26|346.666667|365|0.063142|0.641406|0.000001|
|fallback_only|83.666667|117|0.061133|0.632605|0.001224|
|anchor_only|228.666667|282|0.085583|0.767328|0.000001|
|calibrated_anchor|93.000000|137|0.081632|0.733515|0.000001|
|relative|228.666667|282|0.085583|0.767328|0.000001|
|calibrated_relative|93.000000|137|0.081632|0.733515|0.000001|
|panel_relative|93.000000|137|0.081632|0.733515|0.000001|
|target_shuffled|93.000000|137|0.081632|0.733515|0.000001|

Three control splits are cell-disjoint within each case but reuse biology across seeds. NTC fitting cells are separate from templates, actual and reference populations. Guide/batch overlap is archived; this is not guide- or batch-disjoint validation. Null relative descriptors are zero, so relative models share their corresponding anchor populations. Passing this check alone would not establish target-specific treatment accuracy.

## H1 held-family treatment diagnostic

|candidate|pds|fidelity|reach|jaccard|nmae|common_expression_numerator|shape_wasserstein|
|---|---|---|---|---|---|---|---|
|v24|1.000000|0.646629|0.048033|0.057413|0.880268|0.004126|0.382605|
|v26|1.000000|0.593623|0.065601|0.036305|0.880268|0.003815|0.280836|
|fallback_only|1.000000|0.570421|0.041439|0.037764|0.880220|0.003815|0.244253|
|anchor_only|1.000000|0.642887|0.113632|0.034653|0.880268|0.003815|0.318775|
|calibrated_anchor|1.000000|0.636355|0.100806|0.035510|0.880268|0.003815|0.322259|
|relative|1.000000|0.646817|0.115789|0.034305|0.880268|0.003815|0.322863|
|calibrated_relative|1.000000|0.646369|0.106160|0.034974|0.880261|0.003815|0.305687|
|panel_relative|1.000000|0.644011|0.101782|0.035802|0.880268|0.003815|0.308481|
|target_shuffled|1.000000|0.640701|0.104763|0.035172|0.880268|0.003815|0.322790|

These are raw public-development metrics, not anchor-scaled leaderboard estimates. Treatment means remain fixed to the family-excluded V24 algorithm, not the published ABC populations. All query treatment cells are excluded from fitting and inference. In the H1 pilot, all nine arms are evaluated even if null calibration fails; this diagnostic exception was prespecified before outcomes.

Null checks passed: **0**. Failed checks: H1_260260_DE_calls, H1_260260_detection, H1_260260_variance, H1_270271_DE_calls, H1_270271_detection, H1_270271_variance, H1_270272_DE_calls, H1_270272_detection, H1_270272_variance.

Treatment checks passed: **0**. Failed checks: jaccard, three_metrics_improve, beats_calibrated_relative.

## Fitting and audit

|model|epochs|best_validation_ELBO|source_contexts|
|---|---|---|---|
|H1_panel|12|8.396646|HCT116, HEK293T, CD4_Rest_D2, CD4_Stim48_D2|

The panel network permits measured panel outputs while masking the directly perturbed output per example. Stable state, posterior features and response basis remain nonpanel. Same source sampling, architecture, initialization and held-cell ELBO selection as V26; better likelihood is not a transfer score.

Independent audit passed: **1**. Checked 4 source contexts, 128 complete original rows, 662 fit-only group means and 11 paired cases. 2 full-axis final-emitter replays matched count digests and scorer call counts. H1 source arrays match the previously audited immutable V26 training arrays.

## Release decision

Locked finalist failed prespecified control/treatment checks. All nine pilot arms were intervention-tested; no adaptive substitution. Preserve V24.

No new package, upload or reserve opening occurred. Broader-family fitting was not executed after the H1 rejection. No diagnostic arm was promoted after the fixed finalist failed. The old champion stays available; uploading it again would duplicate an existing result.

## Post-freeze projection diagnostic

|stage|DE_calls|detection_MAE|log_variance_MAE|mean_reference_L1|
|---|---|---|---|---|
|global_integer_repair|0|0.027056|0.285920|0.003188|
|global_CPM_without_strata|2|0.028427|0.316904|0.000001|
|depth_strata_repair|239|0.084182|0.734976|0.001315|
|frozen_final_CPM|240|0.084629|0.742405|0.000001|

This diagnostic uses the already-spent first H1 null. Removing depth strata and retaining global CPM swaps meets all first-case envelopes, but is not an eligible V27 arm and does not revise its rejection. It motivates the separately frozen V28 experiment.

## Operational parallelism

After null comparisons completed, only the verified original orchestration process was stopped. Completed null tables/caches/decisions were preserved. Root run_v27_with_processes.py resumes the identical frozen pipeline using four independent joblib processes and two BLAS threads each. test_v27_processes.py verifies serial/process Atlas arithmetic; the full-axis null/treatment replay audit verifies final scientific counts. analysis/orchestration_handoff.json and process_execution.json record the change. No model, threshold, seed, score or finalist changed.

## Limits and reproduction

One held-out source context cannot establish broad generalization. Public development biology has been reused in prior experiments; three NTC seeds do not add independent lineages. Quantile calibration may alter joint covariance and depth dependence. Source states are observational proxies, and relative Gaussian latent draws are not identified counterfactual cells. Multiplicative relative logits cannot create expression from a zero template entry before quantile mapping or moment repair. The external V25 response basis used allowed source summaries including source held-cell validation observations; validation ELBO is not wholly independent of that preprocessing, while complete H1 exclusion remains enforced. Exact mean constraints can themselves create artifacts. The learning/output-mask and calibration ideas must be judged by matched treatment comparisons, not just control reconstruction.

Methods, inputs, release criteria and commands: [V27_STRATEGY.md](V27_STRATEGY.md). Re-run the pinned scientific runtime; frozen outcomes refuse overwrite and verified caches can resume. Run `audit_v27_experiment.py` then `document_v27_results.py` to verify and regenerate reporting. Large arrays/weights/count caches are regenerable and ignored; root code and small evidence are versioned.
