# V5 model strategy and reproduction

## Source inheritance and isolation

V5 is written to `experiments/v5_atlas_calibrated/`. Its parent is the frozen V4 model and source-component arrays, whose hashes are checked. Original competition control hashes match V4. Earlier packages are preserved. The upstream numerical solver is imported from the reviewed AtlasShift clone after its SHA-256 is checked against the pinned review manifest. The clone retains its MIT license. There is no dependency upgrade or execution of the upstream full-atlas preparation pipeline.

The reviewed repositories, commits, initial control reconstruction experiment and limitations are documented in [UPSTREAM_REPOSITORY_REVIEW.md](UPSTREAM_REPOSITORY_REVIEW.md). This candidate implements the supported ideas with explicit source-quality filters and checks.

## Bounded public H1 extraction

`extract_v5_h1.py` verifies the existing 2025 training subset, reads its original cell identifiers, and extracts target/control rows from generation-pinned public validation/test URLs. Its default limit is 4 GiB of used HTTP ranges per added split, including cached blocks. It chooses the largest predefined cells-per-target/batch limit fitting the planned ranges, from 16, 8, 4, 2 or 1. Rows are selected in source order; this is a bounded sampling approximation, not the complete atlas.

All 1,536 selected training controls are matched by identifier in both extra splits and compared across the entire gene axis. Their count equality is asserted before keeping one copy in the combined training set. Combined cell IDs must be unique. Selected response observations have finite, nonnegative integer raw counts and positive total counts. The measured H1 gene axis is retained separately; unmeasured official genes are supplied by inherited response sources, not imputed source measurements.

`source_extraction.json` and `h1_source_coverage.csv` retain source URL generations/ETags, counts, byte limits and original extraction hashes. Cached HTTP ranges have SHA sidecars. Source arrays are in ignored `reference/v5_h1/`.

## Source contrasts and target/batch validation

`fit_v5_submission.py` uses the same globally sorted batch partition for every public split. Left/right calibration controls are therefore disjoint by cell identifier. It forms both partition-specific and full-source target contrasts using controls matched to represented batches.

The mean-CPM contrast is `log2((target_mean_CPM + 1)/(matched_control_mean_CPM + 1))`; its SE includes target and control mean variances. Pooled probabilities are sums of raw gene counts divided by summed libraries. For their ratio uncertainty, retain covariance through the residual `count_gene - probability_gene * cell_library`; compute variance of this residual, divide by sample size and squared mean library. The pooled contrast is `log1p(50000 * target_probability) - log1p(50000 * matched_control_probability)`. Control batch contributions are weighted by target library mass. Both contrast spaces are clipped to +/-3. Source truth significant sets use actual two-sided Mann-Whitney/Wilcoxon tests of cell CPM with within-target BH; selection calls use the stated 400-cell normal approximation.

The reused V4 H1 selector has three outer target folds/two inner folds and selects specific/prior coefficients while guarding zero-relative response MSE at 1.02. Its selected mean-CPM coefficients are `{"prior_gain": 0.0, "specific_gain": 1.5}`. Pooled gain candidates are 0.3, 0.5, 0.7 and 1.0, with five outer target folds. Gain selection must not exceed the training error of the conservative 0.7 reference. Adoption of a changed gain additionally requires no nested error increase and no composite-proxy decline. The selected pooled gain is **1.0**. This conservative reference is not a historical separate-bulk V4 model.

The mean-response choice has a measurable tradeoff. Against the conservative gain-0.7/no-prior comparison on these same public-source folds, nested response MSE is 0.076939 versus 0.068550; direction fidelity proxy is 0.742 versus 0.520, and Jaccard proxy is 0.185 versus 0.162. Selection favors the combined direction/overlap proxy while retaining the zero-relative MSE guard. This is a source-domain comparison against a conservative estimator, not a hidden V4/V5 score comparison.

Direct-response QC and the other-target shared H1 prior use only qualified targets. Each target is excluded from its own shared prior. The H1 lineage selector operates on source contrasts before the final strict direct-target QC; its diagnostics are source-domain proxies and do not prove anonymous-context accuracy. Complete folds, IDs and candidate scores are saved.

## Context policies and promoter calibration

A keeps V4's matched-lineage CD4 response with unsupported-coordinate K562 fallback. Its pooled effect is the corresponding change in the context's pooled log profile. CD4 gain/centering settings remain `{"prior_gain": 0.0, "specific_gain": 1.0}`; no across-target subtraction is added.

B/C replace supported response coordinates with measured H1 mean and pooled effects. Other coordinates retain `{"cd4_coefficient": 0.5, "k562_coefficient": 1.0, "prior_gain": 0.5}` CD4/K562 transfer coefficients and an other-target H1 shared prior. Pooled transfer from CD4/K562 uses the H1-calibrated gain as an explicit cross-lineage assumption. Source effects are moderated by the unchanged gene and target basal-expression gates `min(1, context_mean_CPM/5)`. Own-target downstream effects are excluded and handled by the separately stated 90% knockdown assumption.

`calibrate_v5_promoters.py` reopens the two disjoint CD4 donor subsets, permits neighboring-KD flags for this narrow calibration, and keeps the remaining on-target/multi-guide/minimum-cell/no-distal/no-low-expression QC. Five outer and three inner target folds select strength from 0, 0.25, 0.5 or 1. The baseline is each source donor's V4-style soft CD4 estimator when eligible, with independent K562 fallback otherwise. The donor-specific prior uses only other training targets. It is applied only if that source donor has significant negative effects in at least two represented conditions; held-out truth-donor results are not used to set a prediction's support gate.

Adoption requires at least ten observed targets, at least 10% nested MSE improvement and at least 70% both-donor negative signs. The accepted fitted close-neighbor log2FC is **-1.663234** (remaining fraction approximately 0.316 before profile normalization). A logarithmic distance ramp is constant through 500 bp and reaches zero correction at 5 kb. Application in A requires repeated significant suppression in full qualified CD4 source conditions and retains stronger existing negative effects. It changed 22 coordinates. Gene-boundary TSS positions are approximate promoter positions, not guide locations; this is not a universal suppression rule or a B/C correction.

## Integer count generation

`build_v5_submission.py` uses seed 20261005. Each context/target has a separate deterministic seed, 400 unique original control rows and balanced source guides/depth quantiles. Templates are individual cell-normalized profiles, with no pooling. Desired mean probability is proportional to the context mean probability times the predicted fold response. Desired pooled probability is obtained by shifting the context pooled `log1p(50000*p)` profile. The own target is multiplied by 0.1 in both profiles before renormalization.

A depth-dependent feasible moment projection bounds the pooled profile while leaving the requested mean fixed. Projection L1 must be at most 0.03. The reviewed iterative solver is capped at 200 iterations; its continuous fitting tolerance/acceptance checks remain active. Integer emission uses seeded systematic rounding and column repair while retaining every original cell depth. The resulting mean probability L1 must be at most 0.03 and pooled probability L1 against the feasible projection at most 0.002. All 900 groups passed. Maximum integer mean residual: 0.007391; maximum feasible pooled residual: 0.000593.

Generation uses a bounded window of four parallel groups by default, with two BLAS threads per fit, and writes groups in official panel order. `--workers 1` provides serial generation with the same seeds/algorithm. Tests verify parallel deterministic output. CSV records retain every original control ID, row, guide and depth plus all group residuals/projections/detection/variance diagnostics. The `CPM_total_variance` diagnostic column stores the sum of variances of cell expression probabilities (counts divided by library), in probability-squared units; multiply by 10^12 to express it as total CPM variance.

## Reproduction from Command Prompt

Prerequisites: `.competition_env`, original competition inputs, unchanged V4 source arrays/model, original CD4 donor/merged caches and H1 training subset, and the reviewed upstream clone. The exact source/clone hashes are in the parent review and V5 manifests. No model evaluation uses hidden 2026 responses or scoreboard tuning.

```cmd
.\.competition_env\Scripts\python.exe extract_v5_h1.py --budget-mib 4096
.\.competition_env\Scripts\python.exe calibrate_v5_promoters.py
.\.competition_env\Scripts\python.exe fit_v5_submission.py
.\.competition_env\Scripts\python.exe -m unittest test_v5_submission -v > experiments\v5_atlas_calibrated\analysis\test_results.txt 2>&1
.\.competition_env\Scripts\python.exe validate_v5.py --model-only
.\.competition_env\Scripts\python.exe build_v5_submission.py --preflight --workers 4
.\.competition_env\Scripts\python.exe build_v5_submission.py --workers 4
.\.competition_env\Scripts\python.exe validate_v5.py
.\.competition_env\Scripts\python.exe document_v5_submission.py
```

The build refuses to overwrite existing final predictions unless `--force` is explicitly supplied. `--package-only` packages an existing H5AD; it does not submit. Cached fitted contrasts require matching script/input hashes. Extracted source reuse requires checksums and pinned object identity; `extractor_code_versions.json` preserves the original and accelerated extractor producers, which use identical row-selection rules and byte limits. If the fitter is intentionally changed, preserve or regenerate its old contrast cache with its provenance rather than bypassing checks. The official VCC 0.2.2 preparation path validates raw counts and packages the file; the private upstream VCC 0.2.0 adapter is not used.

## Verification and limits

Five focused tests cover library-ratio covariance, distinct profile/knockdown semantics, held-target prior exclusion, deterministic integer/parallel generation and infeasible constraints. Independent validation spot-checks newly extracted rows against checksum-verified offline original HDF5 ranges, checks all duplicated selected control counts, recalculates every source mean/pooled effect and uncertainty, and reassembles all 33,359,400 final mean/pooled effects. It also checks target/batch separation, every generated raw count, all 360,000 exact cell depths and original control records, the completed package hash and twelve complete deterministic groups.

All generated group diagnostics and official preparation output are retained. Count distributions are tied to supplied controls: perturbation-dependent cell-depth changes are not modeled. H1 source sampling, selected target panels and inferred lineage matches limit generalization. Neutral profile reconstruction and source-target validation cannot establish the hidden response score; compare the eventual published V5 result with V4 on the same panel/anchors. [V5_SUBMISSION.md](V5_SUBMISSION.md) records the verified artifact and upload command.
