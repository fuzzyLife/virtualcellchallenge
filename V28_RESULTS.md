# V28 results: global transport without depth strata

**No V28 submission is released.** The sole candidate failed prospective checks; V24 remains best recorded at0.14299785521745675.

Frozen plan SHA256: `5cf61e0d9b0031316fa01178f0646a59247d84b398c562d066c5030f7559fa20`. V27 methods/outcomes remain unchanged. This follow-up was selected from spent control/sampler evidence, not V28 treatment scores.

## Real-control calibration

|context|arm|mean_DE_calls|max_DE_calls|detection_MAE|log_variance_MAE|max_mean_L1|
|---|---|---|---|---|---|---|
|H1|v24|982.333333|1128|0.059941|0.872152|0.000000|
|H1|anchor_only|2.000000|2|0.028796|0.339696|0.000129|
|HCT116|v24|670.000000|1013|0.071724|0.314431|0.000000|
|HCT116|anchor_only|2.333333|3|0.024907|0.096513|0.000229|
|HEK293T|v24|490.000000|940|0.077031|0.377376|0.000000|
|HEK293T|anchor_only|4.000000|4|0.028178|0.099577|0.000132|
|CD4_Rest_D1|v24|407.666667|475|0.064247|0.227108|0.000000|
|CD4_Rest_D1|anchor_only|2.000000|2|0.027504|0.089020|0.000004|
|CD4_Stim48_D2|v24|901.000000|1145|0.074435|0.240585|0.000000|
|CD4_Stim48_D2|anchor_only|2.666667|3|0.024572|0.081276|0.000691|

The null splits are cell-disjoint within cases and use independent actual/reference/template/calibration-fit populations. Three seeds reuse biology and are not independent studies. The first H1 case is spent V26/V27 development biology. Every output retains V24 gene totals and template libraries; final CPM drift is measured. These control diagnostics do not prove treatment or leaderboard accuracy.

H1 null gate passed: **True**. Failed checks: none.

HCT116 null gate passed: **True**. Failed checks: none.

HEK293T null gate passed: **True**. Failed checks: none.

CD4 null gate passed: **False**. Failed checks: mean_contract.

## Independent audit

Passed: **True**. Checked15 paired cases and2 full-axis replays. Real-null checks are recomputed independently; hashes, actual-cell labels, unique/disjoint identities, nonnegative integer axes, exact cell libraries/gene totals and CPM drifts are verified. The first H1 global-anchor replay matches the preexisting diagnostic count digest. H1 neural training is unchanged from the prior128-original-row/662-fitting-mean audit.

## Release decision

Sole global-transport finalist failed frozen checks; no diagnostic substitution or reserve opening. Preserve V24.

No new package, upload or reserve opening occurred. Execution stopped at the first failed frozen required stage; the scope tables list only completed comparisons. No old experiment was overwritten and no diagnostic arm was promoted.

## Reproduction and limitations

See [V28_STRATEGY.md](V28_STRATEGY.md) and [V27_RESULTS.md](V27_RESULTS.md). Root code and small evidence are versioned; raw/model/count caches are regenerable and ignored. Run the pinned scientific Python, then `audit_v28_experiment.py` and `document_v28_results.py`. Frozen numerical plans/outcomes refuse overwrite; verified caches can resume.

The development biology is reused. HCT116 lacks NMAE and CD4 D1 treatment populations have200 actual cells. Guide/batch overlap is explicitly archived. The global anchor transports control variation and fixed predicted means; it does not learn identified target-specific stochastic dynamics. Better null calibration does not guarantee improved treatment prediction, source-to-query transfer or leaderboard rank.
