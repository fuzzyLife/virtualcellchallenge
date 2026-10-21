# After V69: tests, failures and a transparent recovery release

V69 fell from V66's published 0.150882188 to 0.095262845. Its intended mean effects were frozen, but the emitted count profiles were not protected; discrimination and expression accuracy both regressed. This run tested twelve hypothesis families (V70-V81) rather than changing one mixture and uploading immediately. No new predictor met its frozen public requirements.

## What the experiments showed

| Versions | Change | Result |
|---|---|---|
|V70-V71|Batch-matched raw-count teacher; remove identity cap; center or confidence-project response|No material expression gain.|
|V72|Recover small batches by aggregating batch-offset count likelihoods|Source coverage improves; predictive error does not.|
|V73|Null-centered count MAP and separate assay/bulk inversion|No material gain; strongest confidence-filtered correction loses to shuffled evidence.|
|V74|Amplify downstream responses in assay space|Expression error worsens by 2.5%-8.3%.|
|V75|Shrink downstream responses while retaining own knockdown|Broader expression error improves 2.36%, PDS falls 0.04781. Rejected.|
|V76|Dilute the whole normalized response, including own knockdown|Broader expression improves 1.73%, PDS falls 0.02682; generated count PDS/reach/Jaccard also worsen.|
|V77|Remove the shared response across targets|No discrimination gain on 55 metadata-selected cases. Remaining600-pair confirmation not run.|
|V78|Restore control variance at exact V66 gene and cell totals|Small changes give weak gains below thresholds; large changes damage DE precision/overlap.|
|V79|Signed-underdispersion-aware molecule thinning plus innovations|Original prospective criteria fail.|
|V80|Empirical zero/detection hurdle and positive count calibration|Original prospective criteria fail.|
|V81|Opposite-direction residual contraction, with expansion controls|Original prospective criteria fail.|

The new batch-teacher initially dropped many tiny matched batches; V72 fixes that coverage defect and records the larger coverage. That defect was in the new teacher, not an established bug in V66. Correcting it did not improve prediction. V75 also exposed that a three-development/one-confirmation split left singleton PDS panels in H1/HCT116. V76 states its reused all-case PDS scope; V78-V81 reserve two cases/context, leaving at least two for confirmation. Failed/partial jobs are retained and are not counted as completed confirmation.

## Exact means are necessary but not sufficient

V78-V81 retain every realized V66 row library and every gene total, hence bitwise-identical raw profiles and PDS. The unchanged integer projector also limits cell-normalized mean drift. Cell distributions can still change Wilcoxon statistics, DE means and the sampling correction used by expression accuracy. The following are equal-study *raw public* metric differences from paired V66, not hosted scaled scores:

|Arm|FID change|JAC change|REACH change|All development checks|
|---|---:|---:|---:|---|
|V78 control10|-0.000364|-0.000031|-0.008719|False|
|V78 control25|+0.000672|+0.000336|-0.001062|False|
|V78 learned25|+0.000082|-0.000316|+0.001062|False|
|V78 query50|-0.094792|-0.021472|+0.012197|False|
|V78 query100|-0.194676|-0.048186|+0.075164|False|
|V79 nb_query|-0.185957|-0.003498|+0.054272|False|
|V79 turnover50|-0.382035|-0.046184|+0.088983|False|
|V79 turnover80|-0.226393|-0.034513|+0.035279|False|
|V79 source50|-0.374255|-0.045150|+0.070134|False|
|V80 hurdle_coupled|-0.224820|+0.004100|+0.069889|False|
|V80 hurdle_independent|-0.193256|-0.005374|+0.029569|False|
|V80 empirical|-0.298294|+0.004562|+0.119403|False|
|V81 contract25|+0.000191|+0.000685|+0.004711|False|
|V81 contract50|-0.012711|-0.005057|+0.002395|False|

Independent audits replayed **410 completed baseline/candidate/control count comparisons** and **3280 raw tied-rank coordinates**, plus all stored BH arrays, four DE metrics and PDS. These are repeated configurations of reused public cases, not 410 independent biological samples. V76 uses a different primary cohort from V78-V81, so differences across their tables are not comparative estimates. Plans, selections, failures, source exclusions and exact input/code hashes are preserved. The individual model criteria were frozen before their outcomes; none was relaxed to promote a failed candidate. Remaining-target tests and second seeds were conditional on success and were not run for failed development arms. V81's control-only guard was specified but not claimed as passed without execution.

## What is ready to upload

`experiments/v82_v66_recovery/prediction.vcc` is a **byte-identical copy of V66**, our best published model. It is 4,073,666,560 bytes, SHA-256 `be899434d81eac45bb8312077fe789088fe49b9b772971d589dcf5414e8f22f1`. The official container validator passes; equality to the previously fully audited payload protects all 360,000 cells and 900 groups. This is a recovery release after V69, **not a new model, new prediction experiment or improvement beyond V66**. The prior V66 public guards were also imperfect and remain recorded as failed. Its published score is evidence for recovery, not a guarantee of current rank or final-test performance. Re-uploading uses an allowance and is only useful if restoring the current leaderboard entry is desired.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v82_v66_recovery\prediction.vcc -m "V66 recovery after V69 - identical predictions" --wait
```

Nothing was uploaded automatically. All original sealed artifacts were preserved. Do not upload a V70-V81 file: no production file was generated for those rejected models.

## Reproduction and next evidence

`summarize_post69_trials.py` regenerates this report from saved result/audit records. `prepare_v82_recovery.py` refuses to copy the champion unless all twelve pilots have completed and failed release, verifies both copies by SHA-256 and calls the official validator. Run it with the official CLI environment. Run numerical experiments/audits with the separately pinned scientific environment; see `MARGIN_VARIANCE_STRATEGY.md` and the per-version strategies/results for commands. Large arrays remain local and ignored by Git; code and compact provenance are committed locally.

The remaining gap is transferable intervention information and validation transfer, not finding another arbitrary temperature. A useful next test needs genuinely unused response data: for example a new CD4 donor with **all CD4 evidence excluded from fitting**, a fixed target/cohort selection, and a mean/decoder crossed comparison fixed before seeing its treatment values. That donor was not downloaded or tested during this run. Reusing the current cases for further tuning cannot establish that a tiny public improvement will transfer to the challenge. This is a recommendation for new evidence, not a claimed successful model or an implemented future experiment.
