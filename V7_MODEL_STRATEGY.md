# V7 modeling and source integrity

Status: implementation and extraction in progress. A package is upload-ready only after `experiments/v7_multisource_adaptive/analysis/independent_validation.json` records a passing audit. No hidden score is available.

## Data and corrected assumptions

The competition still requires 300 CRISPRi targets, 400 raw-count cells per target and context, A/B/C labels and the exact 18,533-gene axis. There are no local competition treated labels.

V7 keeps the measured K562 and CD4 summaries, the expanded 25-target H1 count cache, and samples panel-target cells plus matched non-targeting controls across the complete public Orion HCT116 and HEK293T release. The converted [SLAF Orion release](https://huggingface.co/datasets/slaf-project/X-Atlas-Orion) is pinned to `598aa5442a4ee9a8f1a383a843ab5ce142904d04`. Original Xaira inventory is pinned to `53a5bc98d49247bcf967500292575c3d3602de31`. Sampling is random within each batch and round-robin across batches, capped at 200 treated cells in each of two disjoint batch partitions per target, with up to 32 controls per batch. It covers the full release's target/batch scope rather than every raw cell.

The full K562 raw probe found the same 8,248 measured genes as its pseudobulk, including 7,681 challenge genes. Downloading its 65.8 GB raw matrix would not expand the expression axis. The earlier plan's broader-K562-axis assumption was corrected.

The converted HCT116 expression table omits a block of 5,000 cell IDs while retaining them in metadata; later global row offsets consequently become incorrect. `audit_v7_source_offsets.py` audits every expression fragment's first/last cell ID, identifies gaps, sums missing expression record counts from metadata, excludes selected cells with missing expression and corrects later offsets. It verifies retained cell starts distributed over the complete axis. Every extracted chunk additionally verifies every cell ID, nonnegative whole counts and the complete source library before projection. Verified chunks are reused only if their rows, offsets and counts remain unchanged. Full findings and exclusions are recorded in the extraction manifest. The independent original-Parquet replay also exposed converted duplicate-name suffixes: five extra entries for SPATA13, GOLGA8M, MKKS and ELFN2 were omitted by a literal symbol match. Extraction now uses the original canonical gene metadata by integer token ID and sums duplicate symbols on the challenge axis; all chunks are regenerated under the corrected mapping hash. `validate_v7_original_counts.py` additionally compares 32 complete projected cell vectors and full libraries per line against the original cached Batch1 Parquet. Missing source genes remain masked and do not become measured zero responses.

The reader uses the documented [Lance resolve cache](https://lance.org/guide/object_store/) for public signed-download URLs. `pylance==13.0.0` is installed separately under `reference/v7_runtime`; the competition environment is not upgraded. Source credentials and signed URLs are not included in records.

## Response estimators

`fit_v7_response.py` matches treated source batches to their own controls, computes per-cell CPM with verified source library sizes and constructs log2 effects and delta-method uncertainty. Own-target effects are excluded from response fitting; other panel genes are included. CD4 conditions are correlated, so combining them does not divide their standard error by the square root of three.

`v7_response_methods.py` compares uniform confidence-weighted pooling, pairwise-valid rank-weighted pooling, a target-specific empirical-Bayes mixture shrinker and a conservatively capped measured-knockdown correction. Fixed family weights are K562=2, CD4=0.5, H1=2, HCT116=1 and HEK293T=1. The uniform mode omits extra control-rank weighting; it is not an equal-family vote. Confidence weights are capped; source disagreement contributes uncertainty. The mixture is fit separately for each target using training-source effects. Amplitudes 0.35, 0.65 and 1.0 are tested. Missing source/target/gene estimates abstain.

There are five whole-lineage outer holdouts: K562, CD4, H1, HCT116 and HEK293T. Inner pooling selection excludes both the outer and inner held-out lineage. The production configuration is selected from all public lineage holdouts. The fast response-screen calls are normal/SE proxies; they are not official leaderboard scores.

`train_v7_transition.py` compares three small supervised response heads: a rank-64 linear head, a 64-unit ReLU head and a 128/64-unit ReLU head. Inputs are perturbation identity, frozen RNA-only Tahoe gene vectors projected on training targets, and measured query-control ranks/masks. Loss excludes unmeasured labels, balances source families and excludes the own target. Architecture and training epoch (20, 40, 80) are chosen using an inner held-out lineage. The evaluated outer lineage's treated labels do not enter fitting or selection. This follows the intervention-transition idea in [Arc STATE](https://github.com/ArcInstitute/state), but is not a reproduction of its foundation model. Observational control-expression reconstruction is not used as causal training truth. `train_v7_ridge.py` separately tests penalties 0.01, 0.1, 1 and 10 with the Tahoe/rank kernel and a three-value amplitude grid. Inner selection excludes both the outer and inner lineage. The residual uses only genes observed across targets within a source context; it excludes all panel outputs, which continue to receive pooled responses. Its 25% blend is an additional generated-cell candidate. Full STATE training remains a separate unfinished experiment.

## Generated-cell comparisons

`evaluate_v7_response.py` compares V4-style source transfer, expanded coverage, selected pooling, reduced simulation noise, 10% moment correction, neural responses, a 25% neural blend and a 25% lower-ridge residual blend. The H1 comparator excludes all H1 responses; the HCT116/HEK293T comparators use the prior sources. This is a source-held comparator, not a replay of the original A/B/C upload.

Evaluation targets are selected by deterministic target-ID hash order with at least 200 distinct treated cells (up to 400 observed per target) and an expressed own target, without selecting on their downstream responses. Controls are split within batch into distinct template and reference pools. Exactly 400 prediction templates and 400 reference cells follow treated batch proportions and never overlap. Two fixed seeds are shared across candidates. The calibration and confirmation target sets are disjoint, with public-source coverage determining panel size. Observed treated sample size is recorded and can be smaller than 400; jackknife corrections account for sampling differences. All candidates use the same observed cells and count-test power.

Tests use actual asymptotic Wilcoxon/BH calls on generated integer counts, with the control CPM >5 filter, own-target exclusion, fidelity, reach, overlap and log-FC error. Profiles also provide unanchored discrimination and corrected expression ratio with jackknife sampling correction. Cell variance and libraries are recorded. These metrics follow the locally reviewed competition formulas but do not reproduce hidden panel anchors. Their utility is an explicitly labeled, equal-family mean of six raw metrics.

The calibration winner must retain fidelity, expression error and dispersion within recorded baseline guards. It is accepted only if the frozen winner improves confirmation utility and passes the same guards; otherwise V4 transfer is retained. Public data have been inspected in prior work: this is development evidence, not a pristine blind benchmark or a guarantee of a hidden score gain.

The neutral 400-cell H1 test produced zero false DE calls for all three count generators. V4 and reduced-noise emission preserved exact input counts at zero requested response; 10% calibration retained 0.99879 of template variance. Eight focused scientific tests pass, including empirical preservation of expected counts under reduced stochastic intervention noise.

## Production and audit

The selected model is fitted on all eligible public sources, then applied to each competition context's own controls. The existing directly measured CD4 effects remain the A primary response; expanded transfer fills unsupported coordinates. This is a prior biological transfer assumption, not an identification of A or validation using A treatment labels. B/C use the selected response model. Residual heads abstain for targets absent from all eligible training contrasts; no unseen-target transfer gain is assumed. Gene/target expression gates and separate 90% own-target knockdown remain explicit assumptions.

`build_v7_submission.py` writes sparse raw integer counts in target/context order using 400 distinct guide/depth-stratified control templates per group. It preserves source cell IDs, guide IDs, libraries and seeds, performs preflight, validates the complete H5AD and calls official `vcc prep`. It never uploads. `validate_v7_submission.py` independently checks hashes, source gaps, lineage and batch splits, calibration/confirmation separation, all 900 groups and 360,000 libraries, exact replay of 12 groups and equality of every compressed package CSR value/index/pointer with the H5AD.

## Reproduce in Windows CMD

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m pip install --target .\reference\v7_runtime pylance==13.0.0
.\.competition_env\Scripts\python.exe prepare_v7_sources.py --metadata-only
.\.competition_env\Scripts\python.exe audit_v7_source_offsets.py
.\.competition_env\Scripts\python.exe prepare_v7_sources.py --workers 2
.\.competition_env\Scripts\python.exe fit_v7_response.py
.\.competition_env\Scripts\python.exe train_v7_transition.py
.\.competition_env\Scripts\python.exe train_v7_ridge.py
.\.competition_env\Scripts\python.exe -m unittest test_v7_response -v
.\.competition_env\Scripts\python.exe run_v7_pipeline.py
```

Existing legacy source caches are prerequisites; their hashes and originating scripts are in the prior submission and embedding-pilot records. Numerical caches and large outputs are ignored by Git; the root scripts and small result/provenance records are tracked. Downloads and extraction are resumable with hash checks. Preserve a completed V7 artifact; a model or generator change should use a new experiment version.
