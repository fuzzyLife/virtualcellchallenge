# V29 results: global transport without depth strata

**No V29 submission is released.** The sole candidate failed prospective checks; V24 remains best recorded at0.14299785521745675.

Frozen plan SHA256: `7f8d9e478e3afc35820529efad5eee188fa73b716426567a1c3af199903d6882`. V27 methods/outcomes remain unchanged. The sole candidate is unchanged from V28; the numerical solver was selected from spent control/margin evidence before V29 treatment scores.

## Real-control calibration

|context|arm|mean_DE_calls|max_DE_calls|detection_MAE|log_variance_MAE|max_mean_L1|
|---|---|---|---|---|---|---|
|H1|v24|982.333333|1128|0.059941|0.872152|0.000000|
|H1|anchor_only|0.000000|0|0.027811|0.308536|0.000104|
|HCT116|v24|670.000000|1013|0.071724|0.314431|0.000000|
|HCT116|anchor_only|0.000000|0|0.024464|0.093867|0.000100|
|HEK293T|v24|490.000000|940|0.077031|0.377376|0.000000|
|HEK293T|anchor_only|0.000000|0|0.027685|0.096412|0.000496|
|CD4_Rest_D1|v24|407.666667|475|0.064247|0.227108|0.000000|
|CD4_Rest_D1|anchor_only|0.000000|0|0.026957|0.086486|0.000130|
|CD4_Stim48_D2|v24|901.000000|1145|0.074435|0.240585|0.000000|
|CD4_Stim48_D2|anchor_only|0.333333|1|0.024145|0.078539|0.000495|

The null splits are cell-disjoint within cases and use independent actual/reference/template/calibration-fit populations. Three seeds reuse biology and are not independent studies. The first H1 case is spent V26/V27 development biology. Every output retains V24 gene totals and template libraries; final CPM drift is measured. These control diagnostics do not prove treatment or leaderboard accuracy.

H1 null gate passed: **True**. Failed checks: none.

HCT116 null gate passed: **True**. Failed checks: none.

HEK293T null gate passed: **True**. Failed checks: none.

CD4 null gate passed: **True**. Failed checks: none.

## Paired held-family treatments

|context|metric|V24|global_anchor|difference|
|---|---|---|---|---|
|H1|pds|1.000000|1.000000|0.000000|
|H1|fidelity|0.646629|0.101921|-0.544708|
|H1|reach|0.048033|0.166379|0.118345|
|H1|jaccard|0.057413|0.024805|-0.032608|
|H1|nmae|0.880268|0.880292|0.000024|
|H1|common_expression_numerator|0.004126|0.003815|-0.000310|
|H1|shape_wasserstein|0.382605|0.376452|-0.006153|

These are raw public-development metrics, not host anchor-scaled scores. PDS uses pooled expression profiles and is mathematically protected by the exact V24 gene totals; NMAE is nearly fixed by the CPM contract. Distribution changes mainly affect significance, direction/coverage ordering and jackknife expression-noise correction. Tiny mean-rounding differences should not be interpreted as substantive mean-predictor improvements.

## H1 neural/calibration diagnostics

|candidate|fidelity|reach|jaccard|nmae|common_expression_numerator|shape_wasserstein|
|---|---|---|---|---|---|---|
|v24|0.646629|0.048033|0.057413|0.880268|0.004126|0.382605|
|anchor_only|0.101921|0.166379|0.024805|0.880292|0.003815|0.376452|
|calibrated_anchor|0.101591|0.155302|0.024378|0.880193|0.003815|0.381240|
|relative|0.102323|0.158763|0.024981|0.880311|0.003815|0.384721|
|calibrated_relative|0.099723|0.147244|0.024335|0.880233|0.003815|0.381503|
|panel_relative|0.101171|0.155185|0.024259|0.880243|0.003815|0.376779|
|target_shuffled|0.100745|0.158697|0.024360|0.880214|0.003815|0.374393|

The global raw anchor is the sole eligible candidate. The relative, panel, calibrated and shuffled arms are matched H1 diagnostics only; they cannot replace a failed finalist. All models reuse the unchanged H1-excluded masked/panel fits.

H1 intermediate treatment guards passed: **False**. Failed checks: fidelity, jaccard, family_fidelity.

## Independent audit

Passed: **True**. Checked23 paired cases and3 full-axis replays. Real-null checks are recomputed independently; hashes, actual-cell labels, unique/disjoint identities, nonnegative integer axes, exact cell libraries/gene totals and CPM drifts are verified. A full H1 null count replay matches the archived Python reference exactly. The stimulated-CD4 probe also matches its preexisting full-axis count digest. H1 neural training is unchanged from the prior128-original-row/662-fitting-mean audit.

## Release decision

Sole global-transport finalist failed frozen checks; no diagnostic substitution or reserve opening. Preserve V24.

No new package, upload or reserve opening occurred. Execution stopped at the first failed frozen required stage; the scope tables list only completed comparisons. No old experiment was overwritten and no diagnostic arm was promoted.

## Retrospective mean-response gap

|target|median_abs_true_log2FC|median_abs_predicted_log2FC|predicted_over_true_absolute_effect_sum|sign_agreement_on_true_DE|raw_NMAE|
|---|---|---|---|---|---|
|HDAC8|0.200398|0.036167|0.218372|0.754820|0.842600|
|NFE2L1|0.173204|0.027648|0.212688|0.614840|0.913034|
|PLAGL2|0.218996|0.037770|0.191779|0.676460|0.872134|
|ACLY|0.192885|0.045848|0.298518|0.662117|0.893303|

On the genuinely DE coordinates of these eight spent H1 cases, predicted absolute effect sums are19--30% of measured effect sums, and raw NMAE is0.83--0.92 (a zero response gives1). This is conditional on observed DE selection and is not an unbiased estimate for every competition gene/context. Retrospective oracle scaling is recorded separately: it helps some targets and worsens ACLY, so a universal multiplier is not justified. Oracle coefficients use observed treatment values and are NEVER deployed.

## Conditional decoder mechanism diagnostic

|target|model|query_condition_norm|median_source_condition_norm|cell_interaction_RMS|cell_interaction_energy_fraction|
|---|---|---|---|---|---|
|HDAC8|masked|0.265496|3.163091|0.001754|0.406694|
|NFE2L1|masked|0.285302|3.163091|0.003767|0.463541|
|PLAGL2|masked|0.522885|3.163091|0.006133|0.537807|
|ACLY|masked|1.259053|3.163091|0.006203|0.232403|
|HDAC8|panel|0.265496|3.163091|0.001722|0.364340|
|NFE2L1|panel|0.285302|3.163091|0.003929|0.460446|
|PLAGL2|panel|0.522885|3.163091|0.006352|0.534113|
|ACLY|panel|1.259053|3.163091|0.006413|0.253924|

Query response-descriptor norms are0.27--1.26 versus a median3.16 for nonzero source fitted-group descriptors. Cell-dependent logit changes have RMS about0.0017--0.0065. The source model receives observed fitting-group means and source response descriptors during training; inference receives legacy predicted means/descriptors, while moment projection protects those legacy means. It therefore does not learn a new query mean-response predictor in this experiment. These diagnostics identify conditioning/mean limitations, not a unique causal explanation or proof every neural architecture fails.

## Batch/depth check

|target|seed|treated_batches|pooled_NTC_batches|median_treated_depth|median_pooled_control_depth|median_log_variance_ratio_actual_pooled_NTC|
|---|---|---|---|---|---|---|
|HDAC8|2026101411|48|48|57617.500000|54876.000000|-0.044995|
|NFE2L1|2026101411|48|48|54816.000000|54876.000000|0.002279|
|PLAGL2|2026101411|48|48|53732.500000|54876.000000|-0.009066|
|ACLY|2026101411|48|48|55449.000000|54876.000000|-0.007707|
|HDAC8|2026101427|48|48|56136.500000|54876.000000|-0.026191|
|NFE2L1|2026101427|48|48|54571.000000|54876.000000|-0.006936|
|PLAGL2|2026101427|48|48|53501.500000|54876.000000|0.003649|
|ACLY|2026101427|48|48|56281.500000|54876.000000|-0.031726|

All four held-out targets span the same48 batch labels as controls. Median log-variance ratios correspond to treated/control ratios about0.956--1.004 over eligible nonpanel genes, with similar library depths. A large overall batch mismatch or global variance collapse is not observed in this pilot. Medians can hide changes at specific DE genes; these retrospective checks do not establish causal variance invariance. Actual treatment metadata is used only in this descriptive analysis, never query inference.

## Next justified test

Train a mean-response head from real source perturbations using the same target/context inputs at fitting and prediction, then couple it to a calibrated distribution model. The current count decoder is supplied true source fitting-group means and cannot improve the frozen V24 query means itself. Effects need better gene selection, direction and strength; averaging additional predictors or multiplying every effect is not supported by these results. Keep this validated count-repair baseline and require whole-family treatment gains before opening the donor reserve.

## Numerical backend verification

Three invariant tests,32 varied compiled/Python parity fixtures and a blocked-pair carrier test passed. The full400x18533 stimulated-CD4 probe matches the original Python count digest exactly. Its initial single-case repair timings are70.7s Python and5.82s compiled. These timing/contract results do not measure biological improvement. NumPy/scientific runtime, Numba, llvmlite and joblib versions are pinned; the kernel disables fast-math and parallel floating reductions.

## Reproduction and limitations

See [V29_STRATEGY.md](V29_STRATEGY.md) and [V27_RESULTS.md](V27_RESULTS.md). Root code and small evidence are versioned; raw/model/count caches are regenerable and ignored. Run the pinned scientific Python, then `audit_v29_experiment.py` and `document_v29_results.py`. Frozen numerical plans/outcomes refuse overwrite; verified caches can resume.

The development biology is reused. HCT116 lacks NMAE and CD4 D1 treatment populations have200 actual cells. Guide/batch overlap is explicitly archived. The global anchor transports control variation and fixed predicted means; it does not learn identified target-specific stochastic dynamics. Better null calibration does not guarantee improved treatment prediction, source-to-query transfer or leaderboard rank.
