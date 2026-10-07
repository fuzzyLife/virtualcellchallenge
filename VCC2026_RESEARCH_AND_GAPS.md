# What stronger VCC2026 entries do, and what we are missing

Research and code checked on 2026-10-07. This review separates public evidence, our measured results, and proposed experiments. It does not infer the undisclosed methods of every higher-ranked team.

## Our low score has identifiable components

V7's published overall is **-0.051386751**; V4 remains best at **-0.050673910**. V7 reused exactly V4's response-effect arrays. Its changed random draws did not test the expanded-source estimator on the leaderboard. Repeatedly emitting the same biological model has therefore produced little progress.

The largest negative component is direction fidelity: scaled **-0.672675846**, contributing **-0.112112641** to V7's six-component average. Fidelity yield depends on both correct directions and response coverage; the aggregate alone cannot distinguish too few calls from too many incorrect calls. Its raw value is **0.312449243**. Raw significant-set Jaccard is **0.019646387** and raw corrected expression ratio is **1.245516113**. Expression's displayed zero is a clamped score, not perfect prediction. Improving fidelity to its baseline alone would add 0.112113 to overall; this is arithmetic, not a prediction of achievable improvement.

Full metric comparisons and immutable publication records: [V7_LEADERBOARD_RESULTS.md](V7_LEADERBOARD_RESULTS.md).

## What the public evidence actually shows

[AtlasShift](https://github.com/kaipengm2/Virtual-Cell-Challenge-2026) reports **0.1545618019**, rank 82 on September 10. It combines five prepared public sources and emits counts matching mean CPM and pooled profiles. Its local code separately transfers log2-fold changes for cell means and log1p pooled-profile differences for bulk expression. These are author-reported historical results, not its verified current rank.

[Sidechain's experiment log](https://saberhq.com/sidechain/models/) reports **0.1626** on October 3. Its changes include pooling before target-wide shrinkage, measured knockdown calibration, bounded source confidence, and separate cell-mean/pooled-profile emission. Switching the pooled baseline to depth-weighted controls raised a reported score from 0.1131 to 0.1353; subsequent shrinkage improved expression accuracy. Its neural calibration entry scored -0.9807, and Tahoe-neighbor smoothing barely improved performance. These reports support controlled calibration experiments, not an assumption that larger networks will win. They are participant reports, not independent reproductions.

[Arc's October 1 data article](https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge-2026) explains that 2026 tests unseen cellular contexts using only unperturbed input cells, unlike the 2025 few-shot task. It describes consistent Flex profiling and selection for at least 80% median target reduction. Public screens can differ in platform, timing and knockdown efficiency. Identifying a similar lineage from control markers does not establish a matching intervention response.

The [VCBench preprint](https://arxiv.org/abs/2604.27646) finds that generalization degrades under genuine context shifts, naive dataset merging can hurt, and methods rank differently across metrics. This broader benchmark is not a VCC2026 leaderboard result. [Arc STATE](https://github.com/ArcInstitute/state) provides a genetic-perturbation transition-model example and whole-context holdouts; our control reconstruction networks and small response heads are not a reproduction of that model.

Retrieval hashes, URLs and upstream commits are recorded in [research_sources.json](analysis/transfer_validation_audit/research_sources.json). Sidechain: `da311d9348d6c6f6b33c206bf48faa6c334346ac`; AtlasShift: `d24ce4fdae0cd7cba3ba29546cd8737094eae98a`; official cell-eval2: `0ce25cf4c495e76e60dfca4ddc8f822507aa0e18`.

## A proven problem in our validation

The shared `v6_methods.summarize_profiles` receives one matched control profile per target. Subtracting these different profiles produces artificial target-specific vectors, even when every prediction is identical. The [official metric definition](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics.md) uses one common reference control profile and requires identical predictions to score exactly 0.5.

The new root script [audit_transfer_validation.py](audit_transfer_validation.py) verifies all 94 saved V7 profile-cache hashes, reconstructs disjoint control populations, and compares both origins. Identical predictions score **0.854545 to 1.000000** under historical local semantics and **0.500000** with a common control. This is a demonstrated evaluation flaw, not evidence that our model predicts perturbations perfectly.

The audit also identifies a separate design limitation: V7 prediction templates were sampled using the actual treated cells' batch proportions. Those proportions are unavailable for challenge targets. Even a corrected common-control score can reflect batch composition, especially with small public panels. Consequently, this retrospective audit is not a replacement model-selection benchmark.

Historical numerical results and submitted packages remain intact. Expanded-source candidates also failed the expression guard; the PDS flaw alone does not justify deploying them. See [TRANSFER_VALIDATION_AUDIT.md](TRANSFER_VALIDATION_AUDIT.md).

## Already tried versus genuinely missing

| Area | What we actually did | Remaining gap |
|---|---|---|
| Genetic response coverage | V7 extracted HCT116/HEK293T panel targets across the releases; retained 282/291 eligible contrasts | Expanded estimates were rejected locally and never became the published V7 model; selected-cell summaries are not full-data statistics |
| Output gene coverage | K562 overlaps 7,681 of 18,533 challenge genes; the inspected full raw GWPS axis does not expand that | Many responses rely on other sources or conservative fallback; absent measurement must retain an explicit mask |
| Pooling and shrinkage | V7 tested bounded confidence, source disagreement, adaptive pooled shrinkage and capped knockdown correction | Adaptive shrinkage was tested with rank-weighted pooling, not independently on the selected uniform pool; uncertainty assumptions need calibration |
| Pooled expression | V5 tested dual moments while changing several other components; V7 tested a 10% moment correction | A fixed-response, full decoder comparison separating mean-CPM and bulk effect spaces has not been isolated under a trustworthy benchmark |
| Dispersion | V7 reduced added intervention noise while retaining control-template heterogeneity | This is different from calibrating the entire emitted cell distribution; assess detection, variance and DE calls together |
| Neural models | Control reconstruction, intervention ridge, small conditioned response heads, and blends were tested | A genetic-perturbation transition model with proper context holdouts is a separate experiment; good observational reconstruction is insufficient |
| Evaluation | Nested whole-lineage training folds and separate development targets were implemented | Profile scoring used varying control origins; batch matching used target-specific held truth metadata; local confirmation was already inspected |

## The next experiment, in order

1. **Rebuild validation before optimizing.** Define context/batch strata from control and experimental metadata without using observed response values. Hold out entire source lineages. Use one shared reference-control population per evaluation context and an independent template population, with target-independent sampling weights. Restrict comparisons to compatible technical strata or report batch-stratified results. Keep genuine cell counts where 400 independent cells are unavailable; do not manufacture sample size by treating duplicates as independent. Freeze new calibration/confirmation assignments before fitting. Previously inspected targets are development data.
2. **Require null and replicate checks.** Identical profiles must have PDS 0.5; independent control resampling should show chance discrimination within sampling variation. Record same-context held replicate recovery, control-only expression error, denominator stability, false DE calls and the measured effect size. Use the pinned official scorer for count-level checks and the full common-control gene filter/BH scope for DE.
3. **Run an isolated transfer factorial.** Keep source selection and decoder fixed. Compare uniform pooling with and without target-adaptive shrinkage, then add measured knockdown normalization separately. Include shuffled source-relevance weights as a negative control for learned/context-based weighting. Tune on public calibration only and confirm across lineages and paired seeds.
4. **Run an isolated decoder factorial.** Keep the winning response fixed. Compare the current emitter with pooled-baseline correction and full mean/bulk matching, first at identical response strength. Subsequently test separate bulk strength and total dispersion. Measure requested versus realized profiles and all six metrics; a PDS gain must not conceal expression or DE degradation.
5. **Add neural response residuals only after a strong baseline.** Compare with the same sources, decoder and held contexts. Train on intervention labels, not control covariance. A STATE experiment needs genetic training, gene-axis handling and raw-count emission before it can be compared fairly.

A new submission should represent a changed, validated response estimator or decoder. The current evidence supports fixing the evaluation and performing these controlled tests; it does not establish a new model that will beat V4.
