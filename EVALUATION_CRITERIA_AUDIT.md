# Evaluation criteria audit â€” 2026-10-09

We have a real local-scoring omission and a larger validation-design gap. The official scorer revision is unchanged since our V8 review. Our frozen V11 package is valid and its published score is authoritative; local proxy scores should not be treated as forecasts of that score.

## Sources and scope

Read the [live evaluation page](https://virtualcellchallenge.org/evaluation), checked its public evaluation JavaScript against the supplied screenshot, and compared the public scorer with our code. The live page confirms equal weighting of six metrics across three contexts, integer counts, 400 cells per group, 18,533 ordered genes, no non-targeting rows, a 4.75-billion stored-entry cap and a one-million-count cell-depth cap.

Official `main` is `0ce25cf4c495e76e60dfca4ddc8f822507aa0e18` (2026-10-01), the same revision already pinned in `analysis/v8_official_metric_review.json`. All nine previously downloaded source files match byte for byte. Small additional official files and their MIT license are stored under `analysis/evaluation_criteria_audit/official/`; `sources.json` records URLs and hashes. The scorer package identifies version 0.18.0; documentation examples describe reference bundles built using 0.15.0 and rule version 3. These are not evidence of a hosted scoring change. V11's returned anchor ID ends in `r4`; an anchor artifact version is distinct from the rule version. The API record does not identify the service's numerical backend or expose all 18 context-specific anchors.

Metric-layer tests execute the unchanged official DE modules on synthetic DE tables. Namespace-package loading bypasses unrelated package initialization, not metric code. This audit does not execute the hosted service, the full official count-to-DE pipeline, or official hidden reference data.

## 1. Confirmed omission: fold-change gates smaller than 10 genes

Official NMAE excludes the perturbed gene, keeps finite reference-significant fold changes, and omits perturbations with fewer than 10 eligible genes or a zero/nonfinite reference denominator. Missing/nonfinite predicted fold changes count as zero; they cannot shrink the gate. Our `v6_methods.direction_metrics`, imported by the V11 evaluator, scores any nonempty gate and does not guard a zero denominator. A nine-gene fixture returns local NMAE 1 but no official NMAE value; ten genes returns 1 in both. See [official NMAE implementation](https://github.com/ArcInstitute/cell-eval2/blob/0ce25cf4c495e76e60dfca4ddc8f822507aa0e18/src/cell_eval2/metrics/de.py).

Auditing the saved V11 winner (`consensus20_noise65`) gives:

| Stage | Control reference | Cases | Scored locally with 1â€“9 DE genes | Zero DE genes | Eligible with official gate | Eligible source families | Old local NMAE | Gate-10 NMAE |
|---|---|---:|---:|---:|---:|---|---:|---:|
| development | full_independent | 36 | 15 | 15 | 6 | H1, HCT116 | 0.776321 | 0.801975 |
| development | sampled400 | 36 | 14 | 18 | 4 | H1 | 0.754057 | 0.880938 |
| confirmation | full_independent | 28 | 16 | 6 | 6 | H1, HEK293T | 0.893199 | 0.900565 |
| confirmation | sampled400 | 28 | 17 | 6 | 5 | H1, HEK293T | 0.728144 | 0.830252 |

Cases include context/target/seed combinations, not independent biological replicates. These descriptive reaggregations use the existing equal-family/equal-context scheme over available finite values. When a family has no eligible target it drops out: notably, sampled-400 development becomes H1-only. The changed mean therefore changes the measured cohort and cannot be interpreted as a model change or a predicted leaderboard penalty. Neither the original selection nor confirmation records were overwritten. Complete results for all candidates are in `archived_nmae_gate_audit.csv`.

The scientific issue is power: many external perturbations have too few treated cells or too weak a response to support the official fold-change metric. More small-gate tests will not repair that. A future benchmark needs enough independent treated cells and a powered real-side DE cohort in each held-out source family; report supported target counts alongside every metric. Freeze any added targets and their role before comparing candidate models. Keep weak/null perturbations for false-positive diagnostics rather than pretending they supply reliable fold-change accuracy.

## 2. The local selection utility is not the competition objective

The official context score first scales each raw metric `u` using that context's measured mean-response baseline `b` and split-half replicate anchor `r`: `(u-b)/(r-b)`. It then averages the six scaled metrics and the three contexts equally. The expression component is clamped to [0,1]; NMAE has a scaled floor of -6; the other four components are unclamped. A zero expression score means no positive scaled credit, not zero expression error. See [official metric specification](https://github.com/ArcInstitute/cell-eval2/blob/0ce25cf4c495e76e60dfca4ddc8f822507aa0e18/docs/vcc2026_metrics/vcc2026-metrics.md) and the [evaluation page](https://virtualcellchallenge.org/evaluation).

`select_v11_submission.py` instead combines transformed raw metrics with heuristic clipping, balances public source families/strata, and mixes full-independent/sample-400 references at 70/30. This is explicitly an unanchored proxy, and our later research plan also has raw-metric guards. They are useful diagnostics but do not reproduce official relative sensitivity or context baselines. The baseline is an average **perturbation response**, not merely untreated controls. A correct control reconstruction or good baseline gene-expression predictor does not demonstrate perturbation skill.

Future selection should use official public-source anchors generated within each validation fold, plus raw-metric guardrails. Each split-half leg needs its own controls. NMAE has a special anchor: keep the full-reference DE gate and denominator while comparing half-data fold changes; rebuilding the gate independently in each half is not the official estimator. Private A/B/C anchors cannot be reconstructed from released controls alone. Do not substitute the approximate anchor ranges printed in the paper for the actual context-specific values, or scale averaged raw metrics after averaging contexts.

## 3. Fidelity, reach and overlap reward different things

Fidelity is `correct_direction_calls / max(predicted_adjudicable_calls, real_significant_calls)`. The reference need not call a gene significant for a correct predicted sign to count. This is neither ordinary DE precision nor percent-correct over all measured genes. A tested fixture with four real DE genes and 14 directionally correct predicted calls gets fidelity 1 but Jaccard 4/14. We already implement this main formula correctly. See [official direction functions](https://github.com/ArcInstitute/cell-eval2/blob/0ce25cf4c495e76e60dfca4ddc8f822507aa0e18/src/cell_eval2/metrics/direction.py).

Reach ranks reference-significant genes by the prediction's significance, adjusted p-value, raw p-value, absolute fold change and gene name. It takes the deepest prefix with at least 90% correct directions. A leading mistake followed by nine correct genes can recover to full reach. Predicted-nonsignificant genes stay in the ranking pool: the tested four-predicted-call/14-real-response example gets reach 1 when all signs are correct. Therefore increasing DE call counts alone does not assure higher reach. These ranking and purity details are already handled locally.

There is also an uncommon local edge-case mismatch: if reference-significant genes exist but none has an adjudicable sign, official raw reach is 0 while the local code gives NaN. Nonfinite reference LFC handling differs too; current finite-count/epsilon LFC production avoids infinities. These edge cases should be handled explicitly in a new metric layer without changing SHA-pinned historical modules.

V11's published components make our current limitation visible:

| Component | Scaled leaderboard value | Contribution to overall (value/6) |
|---|---:|---:|
| Discrimination | 0.638267 | 0.106378 |
| Expression accuracy | 0.024946 | 0.004158 |
| DE fold-change accuracy | 0.081058 | 0.013510 |
| DE direction fidelity | -0.015859 | -0.002643 |
| DE direction reach | 0.080022 | 0.013337 |
| DE significant-set overlap | 0.013590 | 0.002265 |

The six contributions reproduce overall **0.137003974**. Discrimination contributes **77.6%** of that net total. Fidelity remains below its mean-response baseline; expression and overlap are only slightly above theirs. These are aggregate score contributions, not causal attribution to any modeling component, and the available record cannot locate the weakest A/B/C context.

A +0.2 improvement in any one scaled metric contributes +0.033333 to overall, equally. Improving directionally reliable response coverage and calibrated gene distributions is therefore more promising than optimizing discrimination alone. It is still a hypothesis to test, not a guarantee of a better submission.

## 4. Cell distributions matter even though cells are not paired

The six metrics use pooled expression profiles and per-gene differential-expression tables; they do not require aligning each generated cell to an individual measured cell. Expression/PDS use log1p profiles built from **pooled counts** normalized to 50,000; DE uses arithmetic means and Wilcoxon ranks of **per-cell CPM** normalized to one million. Depth/composition dependence means these are different functionals. V11's dual-moment emitter already addresses their mean discrepancy.

But matching those two means does not match gene-wise dispersion, detection frequency or Wilcoxon ordering. Our prior neutral experiment found a large false-DE burden for the compact emitter even without a biological response. That supports the already-planned controlled distribution experiment: freeze V11 response effects, vary only the emitter, require control-null calibration plus powered held-family perturbation tests and evaluate all six metrics. Gene-gene covariance can improve biological realism, but the scored DE tests are gene-wise; a better control embedding or gene-reconstruction neural network alone is not evidence of a better intervention predictor.

Only control-expressed genes with mean CPM >5 enter the four DE metrics. Strong induction of a previously silent gene can still matter to expression/PDS but receives no DE credit under this filter. BH correction runs within each perturbation on the retained universe, **before** excluding its own target from metric calculations. Our DE filter, normalization, BH scope and own-gene exclusion already follow these principal rules. PDS additionally removes **all 300 panel targets** from its expression comparison, rather than just each row's own gene; this is already handled in V8â€“V11 and tested by the constant-prediction score of 0.5.

## 5. A secondary numerical parity issue needs backend testing

Local DE uses SciPy MWU with continuity correction enabled. The official scanpy path uses a tie-corrected normal rank statistic without that correction; see [scanpy source](https://github.com/scverse/scanpy/blob/main/src/scanpy/tools/_rank_genes_groups.py). A small analytic equal-depth count fixture crosses the 0.05 threshold: corrected p/q **0.051332466**, uncorrected **0.046263820**. This nine-cell example demonstrates a possible decision difference, not its frequency at competition sample sizes. The official preset selects a backend automatically, and scanpy/pdex are not installed here. Hosted backend identity and end-to-end parity remain unverified. Use an isolated, version-pinned official CPU environment and compare real count fixtures before treating our p-values as exact service reproductions.

## Next-test amendment, separate from the frozen research plan

1. Build an isolated official-compatible validation layer: correct NMAE eligibility/omission, handle zero-denominator and reach edge cases, and verify count-to-DE numerical backend parity.
2. Audit validation power before candidate selection; require a reported usable fold-change cohort in every held-out source family. Retain all perturbations for metrics where they are defined and retain nulls as negative controls. Do not impose the NMAE ten-gene gate on the other five metrics.
3. Build fold-specific mean-response and replicate anchors using official estimators. Report six scaled components and their flat average, alongside all raw metrics and support counts. Separate sampling-reference sensitivity from the main objective.
4. Test the dispersion emitter at fixed V11 biology; then independently test sign-confidence calibration and source-data additions. Keep the previously reserved 12 count-confirmation targets untouched. These suggestions do not revise or unfreeze `future_experiment_plan.json`; a prospective implementation must register its amendment before reading new confirmation outcomes.

## Reproduction and verification

Run from the project root in Windows Command Prompt:

```cmd
.\.competition_env\Scripts\python.exe audit_evaluation_criteria.py
.\.competition_env\Scripts\python.exe test_evaluation_criteria.py
.\.competition_env\Scripts\python.exe document_evaluation_criteria.py
```

The collector requires internet access and refuses an unexpected upstream scorer revision. The test runs offline after collection. Seven official-function fixtures pass with expected mismatches explicitly asserted; constant-prediction PDS, the continuity-boundary count fixture and the published-score arithmetic also pass. Tests verify all **18 frozen V11 numerical dependencies** and **17 original generation snapshots**. No model was refitted, no prediction package was regenerated or uploaded, and previous research/confirmation outputs remain frozen. `test_results.json` records runtime versions, source hashes, test outcomes, publication arithmetic and provenance checks.
