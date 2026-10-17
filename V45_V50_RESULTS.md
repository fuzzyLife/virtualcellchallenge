# V45–V50: populations, numerical repair and empirical priors

V24 remains the best published submission: overall0.14299785521745675. All new results are raw public diagnostics, not private anchor-scaled leaderboard scores. Public biology has been repeatedly examined; the additional-count confirmation is not blind cell-line validation.

## What was tested

| Version | Question | Outcome |
|---|---|---|
| V45 | Smaller empirical residual and response-strength gating | Rejected / V24 retained |
| V46 | Ensemble whole cell vectors from compact and more variable generators | Rejected / V24 retained |
| V47 | Amplify V24 directions corroborated by independent KOLF targets | Rejected / V24 retained |
| V48 | Project cell mixtures to exact V24 margins | Rejected / V24 retained |
| V49 | KL projection plus floor/ceil bounded integer rounding | Rejected / V24 retained |
| V50 | Generic output-gene response propensity versus fixed randomized checkpoint | Selected randomized_prior |

## Findings

- V45 passes individual loss guards but its average relative gains are negative.
- V46 improves REACH but drifts finite-population gene totals, reducing PDS/expression; some nulls worsen.
- V47 genuine target corroboration fails. Its fixed target-shuffled control improves three raw metrics and passes ordinary development guards. This does **not** establish useful target identity.
- V48 preserves pooled profiles/PDS exactly, but loses CD4 fidelity and overlap.
- V49 makes rounding safer:30 numerical trials preserve exact margins with every entry floor/ceil of the smooth expectation. The previous greedy repair reaches2.409-count deviation in these examples; the new maximum is0.950. This concerns the projection repair, not a proof of a flaw in the original V24 emitter. Real-cell V49 still loses CD4 fidelity and fails the0.0005 CPM tolerance in HCT116 (best0.00061142). It is rejected with52 partial records retained, not called a completed success.

## V50 interpretation

The negative control suggests a generic amplitude prior may help even when target-matched transfer does not. The averaged mobility prior pools signed KOLF response propensity per output gene, leaving the current intervention out. The randomized prior is the exact previously observed shuffled checkpoint from V47. Its discovery is explicitly outcome-informed; it receives no credit for target-specific biology. Both only strengthen existing nontrivial V24 effects, retain signs and zeros, bound gains1..1.25 and keep the original count generator. z scores are confidence values, never expression counts or LFC.

Selection is frozen before opening additional count cases. Only candidates passing every original six-metric/null guard, mean error bounds, at least two raw improvements and positive oriented average gain qualify. V49/V50 are compared under a frozen cross-version rule, then only one candidate is locked. Additional24 count cases at two seeds cannot be used for reselection. D3 remains unopened. A positive public surrogate is not a predicted VCC overall score.

## V50 development raw metrics

| Arm | PDS | FID | REACH | JAC | NMAE ↓ | Expression ratio ↓ |
|---|---:|---:|---:|---:|---:|---:|
| v24 | 0.791667 | 0.581703 | 0.194389 | 0.070027 | 0.912601 | 37.471936 |
| mobility_prior | 0.791667 | 0.583389 | 0.208418 | 0.070151 | 0.913063 | 37.551120 |
| randomized_prior | 0.791667 | 0.580516 | 0.211868 | 0.070268 | 0.907941 | 37.543811 |

## V50 confirmation raw metrics

| Arm | PDS | FID | REACH | JAC | NMAE ↓ | Expression ratio ↓ |
|---|---:|---:|---:|---:|---:|---:|
| v24 | 0.781250 | 0.586198 | 0.230317 | 0.053056 | 0.912759 | 37.640579 |
| randomized_prior | 0.782738 | 0.586614 | 0.249261 | 0.053273 | 0.914064 | 37.700543 |

Selected arm: **randomized_prior**. All public guards pass: **True**. Private leaderboard score: **unknown**.

Upload file: experiments/v50_response_prior/prediction.vcc (4,073,646,080 bytes). SHA256: 4629cfca53b031527823e18792d51902fde9bd6e19e84513c4f2dcfc7945a4cb.
Official packaging passed; not uploaded. The complete compressed payload is decompressed and every CSR value/index/pointer compared with the staged source. Canonical CLI cell-index renumbering is checked explicitly; gene/context/target axes and row order must match.

## Generation-seed sensitivity

| Seed offset | Oriented relative raw-metric surrogate gain |
|---|---:|
| 0 | +3.340% |
| 380000 | -0.673% |

One confirmation seed improves the surrogate and the other declines. The fixed aggregate acceptance rules pass, but the evidence is seed-sensitive; this is not proof of stable superiority or a significant leaderboard improvement.

## Packaging audit correction

The first official archive matched every CSR value, index and pointer, but the custom metadata reader failed on the nullable-string obs index written by the CLI (values/mask group). The predictions remain unchanged. Original attempt receipts and reader are preserved. The retry uses backed AnnData decoding, tested with the actual nullable/categorical encodings and deliberately corrupted labels/counts.

## Remaining limitations

The original compact generator still creates excessive control-null DE calls compared with real control/control comparisons. V50 leaves these nulls unchanged; FID on nulls is not used as evidence of biological success. No new model solves missing context-specific intervention biology. A small public gain or safe numerical solver does not establish a significant rank improvement. Only the private scorer can measure that.

## Reproduce

Use the pinned Python3.12 scientific environment (requirements-v35-science.txt, plus requirements-v49-rounding.txt for the rounding pilots). Use .competition_env only for the official CLI. Earlier V7/V24/V30/V34/V37 source/cache workflows are prerequisites, each protected by receipts. SCIENCE_PYTHON below is a placeholder for that interpreter; commands are Windows CMD syntax. Freeze once in a fresh experiment directory; existing frozen results are preserved.

~~~bat
SCIENCE_PYTHON test_v49_balanced_rounding.py
SCIENCE_PYTHON test_v49_population.py
SCIENCE_PYTHON run_v49_balanced_population.py freeze
SCIENCE_PYTHON run_v49_balanced_population.py prepare
SCIENCE_PYTHON run_v49_balanced_population.py means
SCIENCE_PYTHON run_v49_balanced_population.py null
SCIENCE_PYTHON run_v49_parallel_counts.py development
REM V49 actual run fails its CPM contract; preserve partial evidence.
SCIENCE_PYTHON record_v49_failure.py

SCIENCE_PYTHON test_v50_response_prior.py
SCIENCE_PYTHON run_v50_response_prior.py freeze
SCIENCE_PYTHON select_v49_v50_final.py freeze
SCIENCE_PYTHON run_v50_response_prior.py prepare
SCIENCE_PYTHON run_v50_response_prior.py means
SCIENCE_PYTHON run_v50_response_prior.py null
SCIENCE_PYTHON run_v50_parallel_counts.py development
SCIENCE_PYTHON run_v50_response_prior.py select
SCIENCE_PYTHON select_v49_v50_final.py select
REM Continue only with the locked qualifying finalist:
SCIENCE_PYTHON assess_v50_release.py lock
SCIENCE_PYTHON run_v50_parallel_counts.py confirmation
SCIENCE_PYTHON assess_v50_release.py assess
SCIENCE_PYTHON audit_v50_results.py
SCIENCE_PYTHON audit_v50_raw_scores.py
REM Generate only if confirmation passes:
SCIENCE_PYTHON build_v50_submission.py fit
SCIENCE_PYTHON audit_v50_inference.py
SCIENCE_PYTHON build_v50_submission.py preflight
SCIENCE_PYTHON build_v50_submission.py generate --workers 4
SCIENCE_PYTHON audit_v50_source_rows.py
.competition_env/Scripts/python.exe package_v50_submission_v2.py
SCIENCE_PYTHON audit_v50_delivery.py
SCIENCE_PYTHON finalize_v50_delivery.py
~~~

Root source scripts, tests and small receipts are committed. Large raw data, arrays, models, provenance CSVs and upload packages remain regenerable and ignored by Git. No automatic upload or Git push.

Primary background: [predictive ensembles](https://proceedings.neurips.cc/paper/2017/hash/9ef2ed4b7fd2c810847ffa5fa85bce38-Abstract.html), [heterogeneous perturbation responses / mixscape](https://www.nature.com/articles/s41588-021-00778-2), and [RNAi/CRISPR consensus limitations](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2003213). These motivate controls and distribution tests; they do not validate our particular candidates.

## Submit manually from Windows CMD

The package is ready; the private score remains unknown. The aggregate public checks pass, with seed-sensitive gains.

~~~bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v50_response_prior\prediction.vcc -m "Empirical response-prior calibration v50" --wait
~~~

The final release manifest binds every proof receipt and the package SHA256. No upload or Git push was performed.
