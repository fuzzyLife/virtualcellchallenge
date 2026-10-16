# V37: complementary component ensembles

## Purpose and boundaries

V24 remains the best published score (0.1429978552). V36 (0.1418454817) improved scaled expression accuracy and overlap but reduced reach. A leaderboard score does not identify an optimal coefficient, gene, context, or cell. This experiment does not assign a different model to each leaderboard metric; one integer matrix must satisfy all six measures.

Existing code reviewed before implementation: V24 teacher and emitter, V32 cached family-excluded baseline, V35 retrieval and target-centering, V36 population-view mean/variance, original raw count scorer and independent audit, original streaming writer and official packaging adapter. Original files remain unchanged.

## Fixed candidates

All transferable effects exclude the directly targeted gene. The established count generator applies 90% own-gene knockdown once. Missing evidence retains V24. Effects are log2 response values, never raw integer counts.

| Candidate | Mean response | Cell population |
|---|---|---|
| mean50 | Half V24, half submitted V36 retrieval contrast | V24 residual10 emission |
| consensus | Retain retrieval-contrast correction only when the population-trained mean agrees in sign; use smaller absolute correction | V24 residual10 emission |
| center_contrast50 | V24 plus half of each fixed centered and retrieval-contrast correction | V24 residual10 emission |
| population50 | Equal mixture of V24 and submitted V36 expected effects | Generate both complete populations and choose a balanced random half of rows from each, independent of expression/labels |
| mean50_noise_mix | Same response as mean50 | Initial weights 0.30 expected + 0.65 Poisson + 0.05 empirical, then unchanged dual-moment count calibration |

The last candidate mixes expected and empirical cell variation; it is not a claim of exact V11 count equivalence. These related models are not independent voters. Consensus is a bounded hypothesis, not proof of biological correctness. Population mixtures preserve each template row and its exact library, but may differ from a single averaged distribution in means and tails.

## Validation and selection

`experiments/v37_component_ensemble/analysis/plan.json` freezes code/model hashes and rules before new outcomes. Existing private V24/V36 outcomes and previous public experiments are known. There is no new blind lineage or private response label.

1. Predict with whole H1, whole Orion (HCT116+HEK293T), or whole CD4 excluded from the corresponding response banks and neural models. Use only query controls as inputs.
2. Reuse all 44 raw mean cases to check response errors. Cache exact source rows and source/model hashes, with raw count subsets for replay.
3. Evaluate zero-response nulls with each emitter, then the existing first four treatment cases in each of five contexts (20 cases, six arms including V24).
4. Select using the frozen guards: expression error no worse, reach no worse, PDS loss at most .005, fidelity loss at most .002, Jaccard loss at most .0005, NMAE at most 1.002 times baseline, mean errors at most 1.002 times baseline, no increased null DE calls, at least three of six raw metrics improved. Average cases within context, contexts within biological study, then studies equally. Undefined context NMAE is omitted; an entirely unresolved study fails.
5. Only if an ensemble qualifies, lock the selection and evaluate the remaining 24 raw treatment cases at two fixed seeds. These count diagnostics were not used for selecting the blend; their mean labels and source biology were previously used, so they are not an untouched validation set. Require the same aggregate guards plus at most .02 loss per study for PDS, fidelity, and reach. Do not retune after confirmation.
6. Independently replay DE-set and population-profile arithmetic. Raw Wilcoxon calculations and reach ordering use the inherited frozen scorer. Three targeted boundary tests cover consensus, identity fallback and population row/library preservation.
7. Promote only a passing candidate. If none qualifies, retain V24; do not create a renamed duplicate package. D3 remains unopened. No automatic upload or Git push.

Raw local expression error and NMAE decrease when improving; all six official scaled leaderboard scores increase when improving. Public diagnostic values are not the competition's anchored scores. One seed at development and two at confirmation do not establish biological replication or statistical significance.

## Reproduce

Use Python 3.12 and `requirements-v35-science.txt`. Original frozen output namespaces cannot be overwritten. Large NPZ arrays, predictions, packages and logs are ignored in Git.

```powershell
$py = 'C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe'
& $py -m unittest test_v37_ensemble -v
& $py run_v37_ensemble.py freeze
& $py run_v37_ensemble.py prepare
& $py run_v37_ensemble.py mean
& $py run_v37_ensemble.py null
& $py run_v37_null_corrected.py
& $py run_v37_ensemble.py development
& $py run_v37_ensemble.py select
# Only when selection is not V24:
& $py run_v37_ensemble.py confirmation
& $py run_v37_ensemble.py confirm
& $py audit_v37_ensemble.py
& $py document_v37_results.py
```

Never rerun `freeze` in the saved namespace. Preparation and count cases resume only from checksum-verified completed caches; selection and confirmation decisions refuse overwrite.

## Null-sampling correction

After the initial null results and mean errors were visible, the new prefix-based null sampling was found to differ from the existing frozen three-way control diagnostic. `run_v37_null_corrected.py` reuses the original V34 fitting/actual/reference split and seed. Original records and `null_primary_results.json` remain unchanged; `null_corrected/amendment.json` records known outcomes and the sampling-only change. The active null receipt points to the corrected arrays. Treatment counts, candidates and thresholds are unchanged. Zero-response mean arms share exactly the same population; the population-mixture arm also selects between identical parents. These identical nulls are stored with deduplication evidence; the noise blend is generated separately.

The original split separates actual controls, reference controls and fitting controls; templates intentionally belong to fitting controls. A first adapter assertion incorrectly demanded template/fitting disjointness and stopped before any corrected count was generated. `run_v37_null_initial.py`, the initial amendment and its failure log are retained. The corrected adapter checks the actual intended boundaries and uses exactly the original rows and seed.

## Corrective population-mixture experiment

Initial H1/HCT116 treatment panels exposed a large mean mismatch after row selection: independently calibrated parents do not imply a calibrated mixed population. A single additional hypothesis is frozen separately in `analysis/recalibrated/plan.json`: preserve the same 50:50 row mixture, then rerun the original dual-moment routine (200 iterations, same seed) with the mean50 target probabilities. No new weight grid or relaxed guard. This is explicitly outcome-known correction, not a prespecified sixth primary arm.

`v37_recalibrated.py` implements emission; `run_v37_recalibrated.py` reuses the frozen treatment cases/scorer and exact mean50 effects. It tests the same three-way nulls and all20 development cases. `select_v37_final.py` locks one qualifying candidate across primary and corrective tests using the same relative-gain score before any remaining24 count cases are opened. If none qualifies, retain V24. Any selected arm still requires two-seed additional-count confirmation.

```powershell
& $py run_v37_recalibrated.py freeze
& $py run_v37_recalibrated.py null
& $py run_v37_recalibrated.py development
& $py run_v37_recalibrated.py select
& $py run_v37_recalibrated.py audit
& $py select_v37_final.py
# Only the selected namespace proceeds to confirmation/confirm.
```

The uncalibrated population50 mean-error table uses its half-logFC target as a proxy; the actual emitted row mixture can deviate from that target. Its full count-level results measure this deviation. The corrective arm restores the specified half-logFC moments after mixing.

## Final targeted combination and early stop

The five primary results and the corrected row-mixture calibration results were known before the final recipe was frozen. The primary mean50 noise blend passed all treatment guards but increased control-null DE calls. Final dual-moment calibration reduced control-null calls in all five contexts for the row mixture, but did not pass its treatment guards. The final separately frozen hypothesis combines the noise blend with final calibration, without row mixing. `noise_recalibrated/plan.json` records this outcome-known motivation, exact formula and unchanged thresholds.

Run `run_v37_noise_recalibrated.py freeze`, then `audit_v37_recalibrated.py --namespace noise_recalibrated --seal`, then `run_v37_noise_recalibrated.py null` and `gate`. If any null call count exceeds V24, treatment scoring is skipped and this arm is rejected. Otherwise run its `development` and `select`. Run `audit_v37_recalibrated.py --namespace noise_recalibrated` afterward for effect aliases, source boundaries, metric arithmetic and baseline parity. `select_v37_final.py` selects among qualifying candidates across all three namespaces before additional-count confirmation. No thresholds or coefficients are tuned to the seventh arm outcomes.

Five boundary tests cover the original blends and exactly-once knockdown/library preservation for both final calibration implementations. Effect aliases are sealed before the corresponding corrective treatment outcomes; baseline sufficient arrays are compared exactly with the primary experiment. These tests do not isolate whether any improvement comes from ensembling or an additional calibration pass; no calibration-only causal claim is made.

## Safe numerical replay

The original namespace is sealed. Use `reproduce_v37.py --directory experiments/replays/v37_check --stage all` in the scientific environment to replay into a fresh directory without modifying original outcomes or README. Individual stages are available with `--stage`; `--help` lists them. Model training is not repeated: this requires the saved V35/V36 checkpoints and public raw input/extraction caches; regenerate those first using the V35/V36 strategy if absent. Stages run in fresh processes to isolate module adapters. The seventh arm stops after nulls if its gate fails. This reuses known biology and results and is reproduction, not independent validation. Additional-count confirmation is omitted from the replay runner because no original arm qualified; it must not be presented as blind follow-up.
