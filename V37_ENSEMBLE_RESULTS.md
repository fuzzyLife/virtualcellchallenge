# V37 component ensemble results

**Preferred model: V24.** None of the five primary blends or either corrective calibrated candidate passed the applicable fixed guards. No new upload package was created; the original V24 package remains preferred. Additional count confirmation was not opened because there was no qualifying candidate.

## Development: unscaled public count metrics

PDS/FID/REACH/JAC increase when favorable; NMAE and expression error decrease. Equal contexts within study, then equal studies. Undefined context NMAE is omitted, never zero; entirely unresolved studies fail.

| Arm | PDS | FID | REACH | JAC | NMAE error | Expression error |
|---|---:|---:|---:|---:|---:|---:|
| v24 | 0.791666667 | 0.581702779 | 0.194389218 | 0.070026786 | 0.912600685 | 37.471936358 |
| mean50 | 0.791666667 | 0.582604890 | 0.181746291 | 0.069627676 | 0.908764860 | 37.357041911 |
| consensus | 0.791666667 | 0.578333350 | 0.202511297 | 0.069427273 | 0.913045269 | 37.332401855 |
| center_contrast50 | 0.763888889 | 0.580491298 | 0.220136988 | 0.069447233 | 0.914304999 | 36.808532604 |
| population50 | 0.750000000 | 0.575145347 | 0.226877507 | 0.067888238 | 0.911985350 | 46.324961687 |
| mean50_noise_mix | 0.791666667 | 0.580635357 | 0.201840280 | 0.070299557 | 0.909811632 | 37.366384640 |
| recalibrated_population50 | 0.777777778 | 0.571796394 | 0.224630164 | 0.068212336 | 0.915333470 | 37.419128035 |
| recalibrated_noise_mix | 0.791666667 | 0.577149702 | 0.224363984 | 0.069643170 | 0.912962251 | 37.422751624 |

## Frozen guard failures

| Arm | Failed checks |
|---|---|
| mean50 | reach |
| consensus | fid, jac, three_metrics_improve |
| center_contrast50 | pds, jac, three_metrics_improve |
| population50 | expression, pds, fid, jac, three_metrics_improve |
| mean50_noise_mix | null_no_worse |
| recalibrated_population50 | pds, fid, jac, nmae, three_metrics_improve |
| recalibrated_noise_mix | fid, three_metrics_improve |

## Corrected control-null DE calls

| Context | Actual independent controls | V24 | Noise mix | Recalibrated row mixture | Recalibrated noise mix |
|---|---:|---:|---:|---:|---:|
| CD4_Rest_D1 | 1 | 516 | 569 | 451 | 494 |
| CD4_Stim48_D2 | 0 | 753 | 828 | 691 | 730 |
| H1 | 0 | 1415 | 1445 | 1339 | 1355 |
| HCT116 | 1 | 780 | 847 | 682 | 701 |
| HEK293T | 0 | 355 | 445 | 292 | 343 |

All groups are unperturbed. Hundreds of synthetic DE calls versus very few actual-control calls remain an unresolved calibration problem; reduced counts do not imply correct treatment recovery. Large displayed expression-error ratios reflect small reference denominators in some contexts and are not official scaled leaderboard scores.

## Actual-cell mean errors

| Arm | Relative squared LFC error change | Relative absolute LFC error change | Relative CPM squared error change |
|---|---:|---:|---:|
| mean50 | -0.0313% | -0.0297% | -0.0907% |
| consensus | -0.0492% | -0.0506% | -0.2154% |
| center_contrast50 | -0.3488% | -0.3119% | -1.2204% |
| population50 | -0.0313% | -0.0297% | -0.0907% |
| mean50_noise_mix | -0.0313% | -0.0297% | -0.0907% |
| recalibrated_population50 | -0.0313% | -0.0297% | -0.0907% |
| recalibrated_noise_mix | -0.0313% | -0.0297% | -0.0907% |

## Evidence and limitations

Five boundary tests passed (including exactly-once knockdown during recalibration). Independently replayed 250 DE-set comparisons and 50 population-profile panels, with source-control/reference/treatment separation, array hashes, effect bounds and own-target masking checked. Generation asserted every row library exactly preserved. Raw Wilcoxon and reach ordering remain inherited from the frozen scorer.

A sixth corrective mixture was frozen separately after initial H1/HCT116 panels exposed broken moment matching after row selection. It applies the original dual-moment calibration after the same balanced mixture. This is outcome-known testing with unchanged thresholds. Exact effect aliases and all25 repeated V24 baseline count diagnostics were independently verified; it is not an isolated proof that model ensembling, rather than recalibration, causes any gain. The uncalibrated row-mixture mean-error row uses its half-logFC target as a proxy; full count metrics measure the emitted population.

A seventh final corrective candidate combines the fixed noise mixture with the same final dual-moment calibration, without row mixing. It was frozen after the other six outcomes were known, with identical thresholds and an early control-null gate. If that gate fails, treatment scoring and confirmation are not opened for this candidate. All earlier results remain unchanged.

Initial prefix-based null sampling was corrected to the original frozen three-way control split after mean and initial null outcomes were known. Original null records and results are preserved; the amendment records this sampling-only change. Treatment counts, blend coefficients and guards were unchanged. Zero-response arms with identical emitters share the same numeric null arrays; the noise blend was generated separately.

The public mean/source biology and prior private V24/V36 scores were already known. Development reuses 20 count cases; the remaining 24 raw cases were reserved for count-level confirmation only if a blend passed. There is no fresh lineage, no new independent biological experiment, and no evidence of statistical significance from these small diagnostic differences. D3 remains unopened. No automatic upload or Git push.

Safe outcome-known replay without overwriting originals: `reproduce_v37.py`. Full formulas, guards and commands: [V37_ENSEMBLE_STRATEGY.md](V37_ENSEMBLE_STRATEGY.md). Saved plan, caches, means, nulls, count diagnostics, selection and independent audit are in `experiments/v37_component_ensemble/analysis/`. Root scripts: `v37_ensemble.py`, `run_v37_ensemble.py`, `test_v37_ensemble.py`, `audit_v37_ensemble.py`, `document_v37_results.py`. Large arrays and count files remain outside Git.
