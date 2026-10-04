# V4 model strategy and reproducibility

V4 uses the verified cluster/lineage evidence to separate response-transfer policies while keeping every predicted cell on its original A/B/C baseline. It targets the earlier deficit of downstream response coverage. It is an experimental candidate; no hidden competition responses were available.

## 1. What the new control analysis changes

The supplied 300-gene archive exactly matches all 16,560,000 original counts. Target-only clustering agrees with A/B/C for 99.9819% of controls; the full 18,533-gene clustering agrees for every cell. The broader CD3D/CD3E/TRAC program supports T-cell-like A, B expresses VIM together with keratins, and C expresses EPCAM/CDH1/keratins. These findings support different source policies, rather than treating contexts as mixtures of named external cell types. No inferred biological fractions or external cell counts are submitted. See [TARGET_LINEAGE_INTERPRETATION.md](TARGET_LINEAGE_INTERPRETATION.md).

The five-reference basal fit is no longer used to blend CD4/H1/K562 responses. Its CD4 condition weights are retained only after the existing 50% pull toward an equal Rest/Stim8hr/Stim48hr consensus. That is a regularized modeling assumption, not an identified culture condition.

## 2. Choose the A response estimator using disjoint CD4 donors

Use the same pinned merged and two disjoint donor-group DE releases and quality filters as V2/V3. Require at least 100 cells, significant on-target knockdown, multiple guides, no low-target-expression flag and no neighboring/distal off-target flags. Exclude each target from its own downstream loss and from shared-prior fitting. Only finite, positive standard errors provide measured gene support.

Compare target-specific gains 0.5/1/1.5 and shared-response gains 0/0.5/1/2. The target-specific estimate is:

```text
specific_gain * observed_log2FC * 0.04 / (0.04 + standard_error^2)
```

Shared priors are condition-specific, inverse-uncertainty-weighted means over other training targets, clipped to +/-0.5 log2FC. They use no scored target's response. Estimates are clipped to +/-3 log2FC.

Five outer folds keep all conditions of a target together; three inner target folds choose settings. Both donor directions are scored. Candidate selection maximizes:

```text
direction_fidelity_proxy + significant_set_Jaccard_proxy
    + 0.5 * (1 - relative_log2FC_squared_error)
```

Candidate selection requires squared error no greater than 1.02 times the zero-response error. This guardrail constrains candidate fitting; it is not a guaranteed held-out error bound. The weights are fixed modeling choices, not an exact reconstruction of competition anchor scaling.

Predicted call sets are an explicit **400-cell normal approximation**: rescale donor standard errors by sqrt(source target cells / 400), derive two-sided normal p-values and apply BH within each row's eligible gene universe. Reference CD4 call sets use the other donor group's published adjusted p-values. Direction fidelity counts correct predicted signs divided by the larger predicted/reference call budget, including measured nonzero reference directions. This tests a source response-coverage proxy, not the competition's Wilcoxon statistic. Source donor subsets are enriched for larger responses.

The final setting is specific gain **1**, shared gain **0**. Most nested folds favor specific gain 1; three additionally select a 0.5 shared gain. The nested selection procedure raises direction-yield proxy from 0.02324 to 0.17701 and mean predicted calls from 2.81 to 31.86, while relative squared error rises from 0.98425 to 1.00436. This is an explicit coverage/error tradeoff, not an improvement in every metric. The final fixed setting uses all candidate outer-fold results, so its selected minimum/maximum is optimistic; the nested result evaluates the selection procedure.

For A, prefer supported CD4 gene/target effects. Effective K562 responses fill genuinely unsupported coordinates. Do not mix H1 into A merely because some basal genes correlate.

## 3. Calibrate H1 direct responses and nonimmune transfer separately

Reuse the pinned public 2025 H1 extraction: 1,236 perturbed cells for 13 overlapping targets and 1,536 non-targeting controls across 48 batches. There is no new download. Recalculate mean CPM, batch-matched baseline means, sampling standard errors and log2 ratios with the existing one-CPM continuity value. That value belongs to the response model; the rank-correlation analyses still use raw values without added counts.

Split source batches into disjoint alternating halves. Actual source reference call sets are two-sided asymptotic Wilcoxon tests on CPM with within-target BH and control mean >5 CPM. Estimator prediction calls use the stated 400-cell normal proxy. Three outer target folds and two inner target folds compare direct gains 0.7/1/1.5 and shared gains 0/0.5/1/2, with one-SE soft thresholding of direct effects.

Validation retains **0.7 direct gain, no shared component** for measured H1 targets, the V3 setting. Nested relative squared error is 0.80299. Direct estimates require at least 40 source cells and at least 30% measured on-target reduction; all 13 pass. Source control CPM <=5 does not supply a measured direct gene effect.

For missing H1 targets, separately evaluate CD4 coefficients 0/0.5/1, K562 coefficients 0/0.5/1 and shared H1 gains 0/0.5/1. For each of the 13 held H1 targets, exclude that target from every H1 prior and fit coefficients using the other targets. CD4/K562 measurements of the same gene intervention remain permitted source-context training information. The held target's own H1 response never enters its prediction. Both disjoint batch partitions provide evaluation. This is a small source-context transfer experiment, not B/C validation.

The final fixed coefficients are **0.5 CD4 target-specific effect + 1.0 K562 effect + 0.5 shared H1 response**. K562 effects already carry the existing 0.3 transfer shrinkage and Poisson-noise soft threshold, so 1.0 here does not mean unshrunk K562 transfer. Coefficients are response-estimator scales, not biological cell fractions; missing sources do not cause arbitrary renormalization to full strength.

The final candidate's fixed-coefficient held-target relative error is 1.01191; the nested selection procedure is 1.02142. Its direction-yield proxy is 0.40430 and mean predicted calls are 323.54 versus a reference budget of 354.27. Candidate fitting obeys the 2% error guardrail, while the small held-out error increase is reported. Acceptance follows the improved nested composite, rather than claiming better squared error. The shared prior is based on 13 selected interventions and may not represent the full challenge panel.

B/C use supported direct H1 estimates when available, otherwise the calibrated transfer. The direct H1 component's shared gain remains zero; the 0.5 shared response is used only for targets without direct H1 coverage. A uses the matched-lineage CD4 policy described above. B and C receive different final effects through their different basal expression and target-expression support; neither is declared to be H1.

## 4. Apply expression compatibility and retain original cell variability

Multiply downstream effects by `min(1, context gene mean CPM / 5)` and `min(1, context target mean CPM / 5)`, then clip to +/-3. These are fixed assumptions that moderate weakly expressed genes and targets, not experimentally fitted thresholds for causal activity. A gene below 5 CPM may still have a real biological effect. The expression check uses the original complete context matrix, never a renamed external expression profile.

Keep the 90% direct-knockdown assumption separately from downstream modeling. The perturbed gene is excluded from scored downstream model coefficients; the simulator supplies its direct factor 0.1.

For each of the 900 target/context groups, independently draw 400 distinct real controls, balanced over all 46 non-targeting guides (8–9 cells each) and stratified by full-library depth. Preserve source-row, guide and count-depth provenance for all 360,000 predicted cells. Unlike V3's common control template, every target has a separately seeded draw. No cells are mixed across A/B/C or artificially assigned to an external lineage.

Use the existing binomial thinning for decreases and depth-scaled zero-aware Poisson additions for increases. Predictions remain finite nonnegative integer counts on the exact 18,533-gene axis. This distribution assumption cannot establish a correct biological response by itself.

## 5. Independent numerical and source verification

`validate_v4.py` verifies source hashes, source/target/gene axes and disjoint outer target folds. It independently recalculates every K562 target response from original pseudobulk rows, every mapped CD4 target/context response from original merged donor estimates, and every H1 direct/shared component from the original selected counts and batch annotations. It reassembles all 16,679,700 context/target/gene effects without calling the fitter's assembly functions.

After generation, it checks all 900 groups, all 360,000 original control-row/guide/depth records, package checksums and official package metadata. Selected groups covering low/high response energy, measured H1 and missing CD4 cases receive independent binomial/Poisson gene-count expectation checks, unchanged-gene exact checks and complete deterministic cell-count replay. The generator additionally verifies every stored count and cell total, then runs official `vcc prep`.

`audit_v4_strategy.py` compares six selected targets per context with V3 using actual generated cells and the provided control Wilcoxon procedure. This measures predicted response activity, not true recall, correct directions or leaderboard performance. Its figure explicitly displays the source error tradeoffs.

## Root code and reproduction in Command Prompt

New root scripts: `fit_lineage_submission.py`, `test_lineage_submission.py`, `validate_v4.py`, `audit_v4_strategy.py`. `build_submission.py` adds an isolated `lineage-response` mode and per-predicted-cell source records while preserving previous defaults. Restore the pinned environment, original controls and existing V3 public source caches first.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m unittest test_submission test_cd4_model test_context_ensemble test_lineage_submission -v
.\.competition_env\Scripts\python.exe fit_lineage_submission.py --select-only
.\.competition_env\Scripts\python.exe fit_lineage_submission.py
.\.competition_env\Scripts\python.exe validate_v4.py --model-only
.\.competition_env\Scripts\python.exe build_submission.py --output-dir experiments/v4_lineage_response --reference-kind lineage-response --upward-mode zero-aware --seed 20261006
.\.competition_env\Scripts\python.exe validate_v4.py --prediction-only
.\.competition_env\Scripts\python.exe audit_v4_strategy.py
```

A completed prediction requires intentional `--force` to regenerate. Fitted effect/source-component caches, H5AD and VCC packages are ignored by Git; code, provenance, validation, source-model-selection tables, original-row records, results and figures are tracked. Prior submission files are kept in their original folders.

Completed package details and upload command are recorded in [V4_SUBMISSION.md](V4_SUBMISSION.md).

## Primary scoring source

The official suite excludes the perturbed gene from all metrics and all 300 panel targets from discrimination; direction fidelity penalizes insufficient response coverage. Its expression metric also accounts for finite-cell sampling. [Official VCC2026 scoring definitions](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics.md). Source-derived estimates and local control diagnostics cannot guarantee a better leaderboard score.
