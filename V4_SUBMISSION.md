# V4: published lineage-guided response submission

The new candidate is generated, officially packaged and independently validated. Published as `w5criCkgkCWj1QOgP08E` with overall score **-0.050674**. Rank was 977 at publication and 975 at the recorded API check; ranks change over time. The package is isolated in `experiments/v4_lineage_response/`; the earlier submission files retain their original locations.

| Check | Completed result |
|---|---|
| Package | `experiments/v4_lineage_response/prediction.vcc` |
| Package bytes | 3,309,240,320 (3.08 GiB) |
| H5AD bytes | 2,475,123,088 |
| Shape | 360,000 cells x 18,533 genes |
| Groups | 900, exactly 400 cells each |
| Stored count entries | 2,085,901,209 |
| Maximum cell total | 52,673 |
| Official CLI preparation | Passed; CLI 0.2.2 |
| Regression tests | 23 passed |
| Original source control rows | All 360,000 verified; 400 unique rows per group; 46 balanced guides |
| Independent count expectation checks | 389,193 across 21 groups, with complete deterministic count replay |
| All final model effects | 16,679,700 independently assembled |
| Generation and packaging time | 28.3 minutes |
| Leaderboard result | Published; overall -0.050674 |

Package SHA-256:

```text
b2cd9da4388d43937a016581b91a84f2c005b2f713fb5f0eae4d5fe1737f015b
```

## Published validation result

The first, V2 and V4 scores use the same validation panel and anchor set. No V3 submission score is recorded, so there is no V3 leaderboard comparison.

| Scaled component | First | V2 | V4 | Change from V2 |
|---|---:|---:|---:|---:|
| overall | -0.276290 | -0.223325 | -0.050674 | +0.172652 |
| pds | 0.019644 | 0.103042 | 0.255302 | +0.152260 |
| mse | 0.000000 | 0.000000 | 0.000000 | +0.000000 |
| nmae | 0.016006 | 0.032352 | 0.067539 | +0.035186 |
| fid | -1.665331 | -1.470413 | -0.683855 | +0.786558 |
| reach | 0.033874 | 0.041380 | 0.086849 | +0.045469 |
| jac | -0.061932 | -0.046314 | -0.029878 | +0.016436 |

Overall improved by **+0.172652** versus V2 and **+0.225616** versus the first entry. All five non-expression scaled components improve. Direction fidelity accounts arithmetically for 75.9% of the overall gain; this is a score contribution, not causal attribution to a model change.

Raw expression error worsens from 1.014294 to 1.247733 (+23.0%), although the scaled expression score stays at zero. Zero is a clamped leaderboard score, not zero prediction error. Raw direction-fidelity yield rises from 0.074675 to 0.309235; it is response-budget-normalized yield rather than sign precision. Raw significant-set Jaccard rises from 0.013545 to 0.019731, leaving substantial response-set mismatch.

V4 is the best of the three published project submissions on this validation panel. Overall remains below the zero-score anchor, with direction fidelity the largest negative component. The package changes lineage policy, effect strength, transfer priors and control sampling together; aggregate scores cannot isolate which change helped, which context remains weak, or whether residual fidelity loss comes from missed responses or wrong directions. The next experiments should balance response coverage with expression error using source validation and controlled ablations.

Records: `leaderboard_result.json`, `published_api_status.json`, `leaderboard_comparison.csv` and `raw_metric_comparison.csv` in the candidate analysis folder. `submission_ready.json` preserves the generation-time snapshot; the separate publication record records completed upload and scoring. Upload took 13 minutes 20 seconds. [Official scoring definitions](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics.md).

Record or refresh the read-only API confirmation and regenerate these tables with the root script `record_v4_result.py`:

```cmd
.\.competition_env\Scripts\python.exe record_v4_result.py
```

The supplied publication rank is preserved separately from the current API rank. This command does not upload, regenerate predictions, or alter the package.

## Strategy

A uses supported CD4 responses with a donor-calibrated target-specific gain of 1 instead of V3's 0.5. K562 fills unsupported coordinates. B/C use measured H1 responses for 13 targets; other targets use 0.5 CD4 target-specific effect + the existing 0.3-shrunk K562 effect + 0.5 shared H1 response. These are response-estimator coefficients, not cell-type fractions. Gene/target expression compatibility moderates weakly expressed coordinates. Every group starts from a separately seeded, guide/depth-stratified draw of 400 original controls.

Nested CD4 selection increases the direction-yield proxy from 0.02324 to 0.17701, with predicted response calls increasing from 2.81 to 31.86 per source contrast. Relative squared error changes from 0.98425 to 1.00436. Public H1 direct validation retains the V3 estimator (relative squared error 0.80299). The cross-lineage selection procedure has direction-yield proxy 0.40430 and relative squared error 1.02142; the final fixed candidate has 1.01191. These are source-domain tradeoffs using a stated 400-cell normal call proxy, not competition scores or proof of improved B/C predictions.

All 284 previously covered targets retain modeled downstream effects in A. B/C have response estimates for all 300 through the shared prior, including targets with no direct source measurements. The H1 prior is based on only 13 selected interventions and may not represent the hidden panel. More predicted response genes do not establish correct signs or a better leaderboard result.

Complete derivation, objective, fold separation, source limitations and reproduction: [V4_MODEL_STRATEGY.md](V4_MODEL_STRATEGY.md).

## Completed upload command (Windows Command Prompt)

The command below already succeeded; rerunning it creates another submission.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v4_lineage_response\prediction.vcc -m "Lineage-guided response coverage v4" --wait
```

## Actual generated-response diagnostic

The table compares six deliberately selected targets in each context against the provided controls using the control-mean >5 CPM filter, two-sided Wilcoxon and within-target BH, excluding the own target. It measures predicted response activity rather than hidden response accuracy.

| context | target_gene | version | downstream_DE_calls | up_calls | down_calls | mean_library | cell_count |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A | AKT2 | v3 | 4 | 0 | 4 | 21034.2600 | 400 |
| A | AKT2 | v4 | 48 | 30 | 18 | 21227.4900 | 400 |
| A | SMARCA5 | v3 | 408 | 194 | 214 | 20839.0725 | 400 |
| A | SMARCA5 | v4 | 2225 | 1252 | 973 | 20634.7825 | 400 |
| A | VIM | v3 | 2 | 1 | 1 | 21004.1775 | 400 |
| A | VIM | v4 | 51 | 30 | 21 | 21117.9100 | 400 |
| A | S100A11 | v3 | 2 | 1 | 1 | 21053.1050 | 400 |
| A | S100A11 | v4 | 0 | 0 | 0 | 21215.6475 | 400 |
| A | FXYD5 | v3 | 0 | 0 | 0 | 20990.2725 | 400 |
| A | FXYD5 | v4 | 15 | 6 | 9 | 21029.4125 | 400 |
| A | KIF21B | v3 | 0 | 0 | 0 | 21024.4525 | 400 |
| A | KIF21B | v4 | 0 | 0 | 0 | 20990.6750 | 400 |
| B | AKT2 | v3 | 21 | 4 | 17 | 19926.5625 | 400 |
| B | AKT2 | v4 | 49 | 11 | 38 | 19860.3525 | 400 |
| B | SMARCA5 | v3 | 612 | 311 | 301 | 19963.9000 | 400 |
| B | SMARCA5 | v4 | 998 | 519 | 479 | 20119.7325 | 400 |
| B | VIM | v3 | 6 | 5 | 1 | 19800.0650 | 400 |
| B | VIM | v4 | 12 | 9 | 3 | 19800.9150 | 400 |
| B | S100A11 | v3 | 1 | 0 | 1 | 19937.9450 | 400 |
| B | S100A11 | v4 | 24 | 11 | 13 | 19830.5075 | 400 |
| B | FXYD5 | v3 | 0 | 0 | 0 | 19877.5975 | 400 |
| B | FXYD5 | v4 | 31 | 21 | 10 | 19702.3625 | 400 |
| B | KIF21B | v3 | 0 | 0 | 0 | 19922.3625 | 400 |
| B | KIF21B | v4 | 0 | 0 | 0 | 19997.4600 | 400 |
| C | AKT2 | v3 | 31 | 9 | 22 | 21121.0575 | 400 |
| C | AKT2 | v4 | 42 | 13 | 29 | 21138.6550 | 400 |
| C | SMARCA5 | v3 | 612 | 358 | 254 | 21243.8925 | 400 |
| C | SMARCA5 | v4 | 1042 | 566 | 476 | 21263.6125 | 400 |
| C | VIM | v3 | 10 | 7 | 3 | 21125.7075 | 400 |
| C | VIM | v4 | 3 | 2 | 1 | 21134.9925 | 400 |
| C | S100A11 | v3 | 9 | 6 | 3 | 21119.4150 | 400 |
| C | S100A11 | v4 | 44 | 26 | 18 | 21059.5175 | 400 |
| C | FXYD5 | v3 | 0 | 0 | 0 | 21097.8250 | 400 |
| C | FXYD5 | v4 | 20 | 16 | 4 | 21010.5350 | 400 |
| C | KIF21B | v3 | 0 | 0 | 0 | 21131.1950 | 400 |
| C | KIF21B | v4 | 0 | 0 | 0 | 21093.7000 | 400 |

![V4 strategy and diagnostics](experiments/v4_lineage_response/analysis/v4_strategy_audit.png)

## Code and records

Root scripts are `fit_lineage_submission.py`, `test_lineage_submission.py`, `validate_v4.py`, `audit_v4_strategy.py`, with generation integration in `build_submission.py`. The audit script also regenerates this results document and the README section after package validation.

Tracked records in the candidate's `analysis/` folder include nested source-selection JSON/CSV tables, source hashes, response coverage, control/gene summaries, all 360,000 source-row records, model and package validation, exact count-replay evidence, official preparation output, run settings, V3/V4 response diagnostics and PNG/PDF figures. Large source-component caches, effect caches, H5AD and VCC files are excluded from Git and rebuilt locally. No new public source download or leaderboard submission occurs in these scripts.
