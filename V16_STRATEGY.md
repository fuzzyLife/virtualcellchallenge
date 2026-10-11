# V16: control-preserving transcript transport and source-confidence calibration

V13 remains the best published project submission. V16 tests the two V15 gaps separately: artificial changes to control variability and unsupported transferred responses. Numerical V13–V15 modules, plans and outcomes are preserved. This workflow does not upload automatically.

## Frozen nine-arm comparison

The arms are the original V13 response, gene-calibrated response, and a per-target L1-amplitude-matched uniform response, each with three generators: original compact counts, balanced transcript updates, and noise-preserving residual10. Development compares 252 populations on 14 targets/two seeds: 224 newly generated populations and 28 hash-verified V15 baseline populations. Actual/template/reference IDs and whole control halves are copied exactly from the V15 plan, avoiding the previous sampling-order replay mistake. Source evaluation shares one template sample across targets within each family/seed; production samples controls per target as before. That difference and repeated public biology limit the inference.

Development uses 400 distinct real treated cells and 400 distinct sampled controls per target/seed in H1, HCT116 and HEK293T. References use the complete disjoint right controls. The nine V15-reserved names remain reserved for one finalist: H1 SLIRP/ZNF32/MED25, HCT116 MYRF/CAT, HEK293T DMTF1/CHD6/PHF19/PTMS. They were not opened in V15. These public biological sources have been inspected previously; this is algorithm verification, not blind biological confirmation.

## Calibrating response evidence

`v16_methods.py` retains the exact V13 response pool and its neutral treatment of unquantifiable knockdown. CD4 conditions remain one source family. Approximate normal tests from family-level effects/uncertainties, followed by per-target BH correction, count independently significant source families agreeing with the predicted direction. These are source-confidence features, not proof of significance in anonymous A/B/C contexts.

Eighteen fixed bins combine same-direction source votes (0, 1, >=2), absolute predicted log2 effect (<0.05, 0.05–<0.2, >=0.2) and source sign agreement (<=0.75, >0.75). For each outer excluded family, inner predictions exclude both that family and the family whose labels fit calibration. Training-label targets omit all nine reserved names; predictors may still use known-target measurements in other permitted sources. Nonpanel output genes fit the calibration, with equal family/target/condition weights.

A bin's regression gain is `(sum(prediction * truth) + ridge) / (sum(prediction^2) + ridge)`, clipped to [0,1]. Ridge is fixed at `0.25 * total prediction energy / 18`, anchoring sparse bins to one. Effects >=0.2 supported by at least two significant families retain gain one. Deployment multiplies V13 by `0.5 + 0.5 * gain`, so it cannot reverse signs, amplify effects, or suppress more than half the baseline. Own-target repression stays in the generator. Uniform controls use one scale per target chosen to match the calibrated response's total absolute log2 amplitude; they are also eligible predictors under the same prospective policy.

## Balanced transcript updates

The new balanced generator starts with each original sampled control cell. With original counts `x`, response factors `f` and library `L`, its conditional expected counts are `m = L * x * f / sum(x * f)`.

- Expected decreases remove original transcripts with independent binomial draws.
- Exactly the number removed is redistributed by a multinomial over positive expected increases.
- Every library remains the original integer total; all counts are nonnegative integers.
- With zero response and no own-target intervention, output is exactly the original control matrix.
- Cell-level zero genes cannot become expressed with this construction. It preserves the cell's support and is unsuitable for modeling de novo gene activation.
- No global population-mean refitting follows. Means differ from V13's fixed mean/bulk constraints because cell states, row normalizers and sampling variability remain. Both mean distances are measured rather than asserted away.

This differs from V12 thinning/immigration, which subsequently applied global dual-moment fitting and could change zero-response controls. Balanced transport is an exploratory molecule-redistribution model, not an optimal-transport solver or a mechanistic CRISPR simulator.

## Noise-preserving residual10

The second generator uses `0.25 * expected + 0.65 * Poisson(expected) + 0.10 * empirical_response_counts`. The empirical component is normalized per original cell after applying response factors. Unlike V15's mixture, the Poisson coefficient remains 0.65 before fitting. The original dual-moment solver still follows, so it can change the final variance; final count diagnostics, not the coefficient identity, determine adoption. This generator does not guarantee exact zero-response identity.

## Tests and release gates

Nine scientific invariants test exact null identity, integer libraries/support/replay, expected molecule redistribution, own-target repression, the Poisson coefficient, held-family label exclusion, reserved-label exclusion, bounded calibration and depth-axis validation.

Control-only evaluation generates 18 populations: compact/balanced/residual10, two seeds and three source families. Each uses 400 control templates, 200 other disjoint actual controls, and full independent right controls. Actual200 is descriptive; null release checks compare generated400 with independent right controls. A new generator must not increase null DE calls by >10%, must improve mean absolute log total-variance distance to the reference, and must keep detection MAE within 10% of compact.

Primary development/verification keeps V15's numerical guard thresholds: PDS/fidelity/reach/Jaccard drops <=0.02/0.01/0.02/0.002; NMAE worsening <=0.01; paired expression numerator relative worsening <=2%; worst-family fidelity/reach drops <=0.03/0.05. At least three available metrics, including a DE metric, must improve. Missing NMAE is not filled with zero and prevents claiming complete six-metric verification. Tiny expression denominators are reported; paired numerators control relative worsening. These are exploratory bounds, not statistical noninferiority.

A calibrated arm must also beat its corresponding amplitude-matched uniform arm on at least two of fidelity/reach/Jaccard. Uniform arms can advance if they pass the ordinary guards. Equal metric-improvement ranks choose among passing candidates. One finalist is frozen before reserved verification; a gene-calibrated finalist also retains its uniform control in that verification. Failed verification cannot select another candidate or retune. No private anchor scaling or expected leaderboard rank is estimated.

## Reproduction in Command Prompt

Preserve frozen plans/results; preparation and assessment refuse to overwrite them. Models and count caches resume only if recorded input/code/byte hashes match.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m unittest test_v16_transport -v
.\.competition_env\Scripts\python.exe prepare_v16_experiment.py
.\.competition_env\Scripts\python.exe run_v16_experiment.py models
.\.competition_env\Scripts\python.exe run_v16_experiment.py null --workers 6
.\.competition_env\Scripts\python.exe run_v16_experiment.py development --workers 6
.\.competition_env\Scripts\python.exe assess_v16_experiment.py development
```

The conditional release driver waits for the started model/null/development stages, assesses development once, audits original identities and baseline replay, then runs reserved verification only for the frozen finalist. It fits, generates, packages and validates a new prediction only after verification passes. It never uploads. Resume the driver with:

```cmd
.\.competition_env\Scripts\python.exe release_v16_submission.py
```

For a passing finalist, the driver runs these steps in order (each is also available independently):

```cmd
.\.competition_env\Scripts\python.exe audit_v16_experiment.py
.\.competition_env\Scripts\python.exe run_v16_experiment.py confirmation --workers 6
.\.competition_env\Scripts\python.exe assess_v16_experiment.py confirmation
.\.competition_env\Scripts\python.exe audit_v16_experiment.py
.\.competition_env\Scripts\python.exe fit_v16_submission.py
.\.competition_env\Scripts\python.exe build_v16_submission.py --preflight --workers 6
.\.competition_env\Scripts\python.exe build_v16_submission.py --workers 6
.\.competition_env\Scripts\python.exe validate_v16_submission.py
.\.competition_env\Scripts\python.exe document_v16_submission.py
```

Do not repeat once-only assessment or fitting commands after their records exist. A failed gate ends release and writes the measured outcome; it does not generate another candidate. Root code and small evidence are versioned; regenerable arrays, counts, packages and full sampling provenance are excluded from Git. Final results and the release decision are recorded in `V16_RESULTS.md` and `V16_SUBMISSION.md`.
