# Post-V11 research and tests

Reviewed and tested on 2026-10-09. V11 remains the published baseline: overall 0.137003974, entry `bmIAGveuEvFiYQTvz0Ma`. This work collects small public inputs, diagnoses remaining failures and freezes future tests. It creates no upload candidate and changes no V11 numerical code, model, package or frozen evidence.

## Main findings

The highest-priority missing opportunities are an independently profiled, probe-based genetic source; gene-specific count-distribution calibration; and calibrated context-transfer uncertainty. More elaborate networks are a later option, after those comparisons. The new neutral tests support the dispersion hypothesis but do not establish improved perturbation or leaderboard performance.

1. The original VIPerturb-seq target manifest contains **300/300 targets**. Only **124** are included in its filtered object, the median available count is **47 cells per target**, and none has 400 cells. This may improve source/platform evidence, not solve validation power.
2. Across six A/B/C control splits, gene-specific gamma-Poisson emission reduces mean null DE calls against 2,000 independent references from **1710.7 to 83.5**. Average total CPM variance moves from 0.233 times real-control variance to 1.051. It still makes up to 291 null calls in one split, so it is not fully calibrated.
3. A 25% empirical-control template blend worsens the null bias despite matching population moments. This tested blend should not replace V11. The participant Sidechain likewise reports that another real-control emission change failed on the board after local improvements; the implementations differ and this is not an independent reproduction of its experiment.
4. Whole-family summary diagnostics find weak target-specific transfer and poor reliability in the strongest heuristic-confidence bins. The current scalar pooling uncertainty is not a calibrated directional posterior. Source masks also leave a median of **4677/18,532 non-target outputs** without a usable response per target. These masks include abundance/QC restrictions, not just unmeasured assay genes.

## Scientific papers and public evidence

The source record in `analysis/post_v11_research/literature_review.json` states exactly which sections or primary excerpts were reviewed. Full publisher pages for some Nature papers were inaccessible; the review uses primary abstracts/method excerpts and author code, and does not claim every supplement was read. VIPerturb-seq and the pinned RNA-fingerprinting source describe preprint work. Arc's experiment article and Sidechain's participant log are listed separately from scientific papers. We cannot infer undisclosed methods of higher-ranked teams.

| Primary source | What it changes in our test design |
| --- | --- |
| [Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation](https://pmc.ncbi.nlm.nih.gov/articles/PMC13271886/) | Shared stress and selection effects can inflate control-referenced similarities. Add generic-perturbation and shuffled-target controls; retain official scoring and add separate target-specific diagnostics. |
| [Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines](https://pmc.ncbi.nlm.nih.gov/articles/PMC12328236/) | Compare new networks with strong linear and mean-effect baselines under the same intervention split. This does not establish that every neural architecture fails at this competition. |
| [Benchmarking algorithms for generalizable single-cell perturbation response prediction](https://www.nature.com/articles/s41592-025-02980-0) | Generalization depends on context and task. Source family and technical platform must be held out explicitly; additional drug or cytokine data are not automatically CRISPRi training labels. |
| [Learning single-cell perturbation responses using neural optimal transport](https://www.nature.com/articles/s41592-023-01969-x) | Unpaired population transport is relevant, but our integer moment fitter is not CellOT. Test distributions and covariance as well as population means; any neural OT model needs intervention training and held-context checks. |
| [Deep generative modeling for single-cell transcriptomics](https://pmc.ncbi.nlm.nih.gov/articles/PMC6289068/) | Gene-specific count dispersion and latent cell states are separable. A gamma-Poisson pilot tests dispersion only; it is not a trained scVI model and does not preserve biological covariance. |
| [Characterizing the molecular regulation of inhibitory immune checkpoints with multi-modal single-cell screens](https://pmc.ncbi.nlm.nih.gov/articles/PMC8011839/) | Guide-bearing cells can include transcriptomic responders and nonresponders. Source efficacy scaling may mix response magnitude with responder fraction. Keep all genuine validation cells; do not assume knockout observations apply quantitatively to CRISPRi. |
| [False discovery rates: a new deal](https://pmc.ncbi.nlm.nih.gov/articles/PMC5379932/) | Estimate directional uncertainty rather than interpreting sign agreement as calibrated confidence. Local false-sign rates require a fitted posterior and validated likelihood; our z proxy is not one. |
| [Flexible statistical methods for estimating and testing effects in genomic studies with multiple conditions](https://pmc.ncbi.nlm.nih.gov/articles/PMC6309609/) | Learn heterogeneous sharing patterns instead of only scalar pooling. Mash does not automatically predict an entirely unmeasured context: that needs an explicit training-only mapping and whole-family validation. |
| [Mapping transcriptional responses to cellular perturbation dictionaries with RNA fingerprinting](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458102/) | Control-derived factors and a rank-one perturbation factor separate baseline state from response strength. Fingerprints live in a normalized residual space and are not directly count fold changes. Query mapping requires observed query responses; challenge controls alone cannot provide those. |
| [Genome-wide single-cell perturbation screens with VIPerturb-seq](https://pmc.ncbi.nlm.nih.gov/articles/PMC12919025/) | Probe-based K562 genetic data may reduce platform mismatch and expand evidence. The manifest gives all 300 targets but few cells per target; QC-filtered membership and actual knockdown need checking. Treat this as K562, not an independent lineage. |
| [Behind the data: Virtual Cell Challenge 2026](https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge-2026) | Probe-based challenge profiling and strong knockdown motivate platform-aware source checks. They do not identify anonymous contexts or establish transferability from controls. |
| [Sidechain public model log, SER-15 and SER-14](https://saberhq.com/sidechain/models/) | A real-control emission experiment improved local tests but sharply reduced the reported leaderboard score. Preserve distribution realism as a diagnostic and require perturbation-level confirmation before replacement. |

## Collected data and method code

The root `collect_post_v11_data.py` performs bounded unauthenticated downloads and verifies hashes when reusing them. The original Zenodo manifest additionally matches the publisher MD5. URLs, retrieval times, sizes, hashes and publisher checksums are in `analysis/post_v11_research/data_manifest.json`.

- **VIPerturb-seq:** original record `10.5281/zenodo.18460279`, manifest and target-level coverage/cell/QC table. The genome-wide experiment is K562; K562, HEK293 and HAP1 also appear in a separate small pilot. The genome-wide source must be grouped with other K562 data when holding out a lineage. The three full bins total about 10.2 GB compressed; their shared NTCs must be deduplicated. We have not downloaded the expression objects, verified their gene/count axes or extracted responses. The supplied Seurat RDS format needs an R/Seurat extraction route; no R runtime is available here. The QC-filtered object alone cannot supply all 300 targets. Do not substitute the VIM-enriched screen for the unbiased genome-wide screen.
- **Reactome:** exact downloaded GMT ZIP, derived pathway membership and target annotation coverage. The contemporaneous version endpoint reports **97**; the ZIP is pinned by bytes/hash rather than assuming its mutable URL guarantees a release. There are **2,868 human pathways**, and **184/300 targets** have annotations. Membership is a prior, not an intervention label. Data are CC0 under the [Reactome license](https://reactome.org/license); attribution is retained here.
- **RNA fingerprinting:** author README, MIT license and `R/fit.R`/`R/utilities.R` pinned at commit `4004c8d04c41cdccff22b982411121cc574d5ed1`. The inspected implementation estimates 30 control PCs, removes baseline structure, then estimates a rank-one perturbation factor and filters it. It has not been installed or run. This is materially different from our existing control reconstruction networks and pooled pseudobulk effects. Mapping an observed treated query to a dictionary cannot directly predict an unobserved challenge response from controls alone.
- **Existing public sources:** original Replogle author inventory was refreshed; the V7 cached source response arrays are hash-verified and reused without changes. No newly verified independent lineage has been added to training.

See `viperturb_target_coverage.csv`, `viperturb_inventory.csv`, `reactome_pathways.csv`, `reactome_membership.csv` and `target_annotation_coverage.csv` in the research directory. VIPerturb data are CC BY 4.0 according to the original record; cite Bradu and colleagues and the dataset DOI when using them. Source code is retained with its MIT license. The data collection contains no challenge treated labels.

## Tests executed: distribution and null bias

`benchmark_post_v11_distributions.py` freezes `neutral_test_plan.json` before execution. Every guide is split into disjoint fit, template and reference populations (50/25/25%). All candidates use the same 400 template depths; references contain 400 or 2,000 distinct cells. No own-target knockdown or response is applied. There are four methods, three contexts, two seeds and two reference sizes: **48 comparisons**. All generated counts are positive stored integers and preserve original library depths exactly.

| Method | Description |
| --- | --- |
| real_resample | Independent real control cells; a benchmark, not moment-fitted |
| v11_noise65 | Exact V11 zero-response emitter and integer dual-moment fitting |
| gene_gamma_poisson | Per-gene latent overdispersion estimated from fit controls after counting-noise correction; independent gamma-Poisson template followed by the same moment solver |
| empirical_blend25 | 25% real-cell composition and 75% compact V11 template, followed by the same moment solver |

The two pilots keep requested mean and pooled compositions fixed. Their initialization changes; this is not a trained scVI or CellOT implementation. Independent gamma-Poisson draws do not retain cross-gene biological covariance. The empirical blend's weak variance and high null DE show that adding real-cell structure before moment fitting does not guarantee preservation afterward.

| method | mean_null_DE_calls | max_null_DE_calls | mean_variance_ratio | mean_detection_MAE |
| --- | --- | --- | --- | --- |
| empirical_blend25 | 3766.166667 | 5021 | 0.229141 | 0.095549 |
| gene_gamma_poisson | 83.500000 | 291 | 1.051468 | 0.020776 |
| real_resample | 0.000000 | 0 | 1.020448 | 0.016875 |
| v11_noise65 | 1710.666667 | 2430 | 0.233169 | 0.081986 |

False/null calls here are genes declared DE by the existing asymptotic Wilcoxon/BH implementation with a shared reference CPM filter greater than 5. Detection MAE is the average absolute gene-detection-fraction error on reference-eligible genes. Total variance matching alone is insufficient: B still has substantial null bias under gamma-Poisson. The per-gene diagnostic CSVs show where variance and detection remain wrong. There is no perturbation-level or official anchored score for these new pilots.

## Tests executed: transfer, uncertainty and systematic signal

`audit_post_v11_transfer.py` evaluates **1462 public source contrasts** from five families/seven contexts, excluding the entire evaluated family from pooled responses and shrinkage. It reproduces the ungated V11 consensus-20 response equation; it does not apply the production context-abundance gate or generate counts. Unsupported contrasts abstain (**5** have no evaluable downstream outputs).

Matched-target, generic training-centroid and shuffled-target cosine comparisons use the same common predictor-supported output gene mask in each context and exclude all 300 panel outputs. Two target-label permutations test identity signal. Separate descriptive centering subtracts evaluation-panel centroids solely while calculating a residual metric; no held response enters fitting. This centered metric complements the official common-control metric and is not a replacement scorer.

| family | lfc_mse_ratio | de_sign_accuracy | heuristic_95_interval_coverage | transfer_cosine | residual_cosine_diagnostic | shuffled_target_cosine | generic_train_centroid_cosine |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CD4 | 1.291753 | 0.646357 | 0.982763 | 0.029244 | 0.029089 | -0.000127 | 0.026539 |
| H1 | 1.643104 | 0.631996 | 0.964766 | 0.057114 | 0.050785 | 0.007374 | 0.045496 |
| HCT116 | 1.027863 | 0.790881 | 0.978446 | 0.032641 | 0.034073 | 0.005480 | 0.106754 |
| HEK293T | 1.213072 | 0.770191 | 0.983580 | 0.035023 | 0.034602 | 0.007004 | 0.067444 |
| K562 | 1.065666 | 0.590113 | 0.978289 | 0.043406 | 0.041725 | 0.001472 | 0.043705 |

These are retrospective, unanchored source-summary diagnostics. Normal/SE source q-values define the reference-DE set; they are not generated-count Wilcoxon calls. The MSE ratios above 1 show that this ungated estimator can be worse than zero effects in fold-change MSE even where some directions transfer. This is not a contradiction of V11's published score: the representation, gating and metric are different. A generic shared effect can have a larger cosine while containing no target-discrimination information.

`uncertainty_bins.csv` separates support, agreement and effect/uncertainty bins. The heuristic uncertainty is pooled transfer variance plus held-source measurement error, not an empirical-Bayes posterior variance. Its normal sign-error proxy and observed sign errors among reference DE genes should therefore be compared descriptively, not treated as a validated false-sign posterior or unconditional calibration curve. The very strong bin has many more errors and much poorer interval coverage than its proxy suggests. Future calibration must be fitted on inner training folds and confirmed independently.

## Small hypothesis ruled down: target knockdown composition

`audit_post_v11_composition.py` tests the exact mean-composition algebra for changing only the target: other-gene log2 fold change is `-log2(1 - knockdown * target_fraction)`. The largest 90%-knockdown effect is **0.007879**, for `VIM` in `B`; changing 80% to 90% changes it by only **0.000878**. This is a real indirect effect even though the own target is excluded from scoring, but its small size makes it a lower priority than the measured distribution and source-transfer gaps. This algebra is not a test of absolute RNA abundance or statistical significance.

## Concrete next experiments, in order

The frozen `future_experiment_plan.json` contains the full designs, negative controls and screening limits. All comparisons retain exact V11 as a baseline. Source additions, response estimation and count emission must be isolated before combining them.

| Test | Hypothesis | Design | Negative control |
| --- | --- | --- | --- |
| D1 | Matched Flex genetic data improve source transfer after controlling for source family. | Extract the 300 targets plus deduplicated NTCs from all VIPerturb bins. Verify raw counts, output gene mask, guide knockdown, cell IDs, batch and library axes. Compare original source pool, Flex replacement within K562, and within-family combined summaries at identical V11 emission. Hold out all K562 datasets together; do not count platforms as independent lineages. | Shuffle target labels within the added source; compare a size-matched downsample of existing K562. |
| D2 | Gene-specific distribution calibration reduces expression error without erasing direction/coverage. | Keep V11 effects fixed; compare exact V11 emission, moment-matched gamma-Poisson, and one 25% gamma/75% compact template. Test null controls first, then actual public perturbations at genuine cell counts. Analyze direction, reach and six-metric panel jointly. | Real control resampling with no response; source target-label permutation; own-target knockdown disabled for nulls. |
| D3 | The strongest response confidence proxies are miscalibrated across contexts. | Fit a bounded transfer-heterogeneity calibration and false-sign reliability curve using inner held-family development folds. Compare scalar empirical Bayes with a multivariate/covariance or pathway prior. Covariances must use training effects only; a novel family has no measured effect and requires an explicit prediction rule. | Permute source identity/uncertainty while preserving marginal weights. Report held-family interval coverage and sign error by abundance/support, not only aggregate MSE. |
| D4 | Cell-state and responder heterogeneity bias source pseudobulk effects. | Compare raw source CPM effects with control-factor-adjusted means and a rank-one perturbation-factor pilot inspired by RNA fingerprinting. Include soft responder mixtures; do not delete weak responders from validation truth. Calibrate factors back to mean/pooled count effects before emission. | Random factors of identical rank; NTC pseudo-perturbations; target-shuffled response factors. Query challenge controls alone cannot identify an unobserved treated fingerprint. |
| D5 | Pathway structure improves residual target-specific response rather than shared stress recovery. | Use pinned Reactome membership for a small residual prior after subtracting training-only shared response. Compare zero effect, training perturbation centroid, V11, and pathway residual. Keep official common-control metrics and report descriptive perturbation-centered diagnostics separately. | Size/degree-matched shuffled pathway sets; exclude each target and its held functional group when testing extrapolation. |

There are 12 reserved, globally disjoint, previously unused count-test targets:

| context | targets | available_cells |
| --- | --- | --- |
| HCT116_controlstratum0 | TXNDC16, PATL1, PBLD | 154, 219, 169 |
| HCT116_controlstratum1 | ZBTB34, MYD88, SYVN1 | 186, 150, 157 |
| HEK293T_controlstratum0 | PCGF6, PCGF1, DMTF1 | 178, 178, 223 |
| HEK293T_controlstratum1 | TEX2, NAIF1, TRIP4 | 161, 158, 170 |

Selection used target names and available cell counts in existing control-defined technical strata, not response values. All previous V7-V11 count-test targets are excluded. Public source summaries have already been inspected, so this is not pristine biological blind validation. H1 remains reused development data. Preserve genuine observed counts (often 150-200), and never duplicate cells to claim 400 independent observations. Use both sampled-400 and full independent references, paired seeds and equal family/technical-stratum weighting.

The initial count-stage is limited to four candidates including V11. Freeze the development winner before opening reserved confirmation. Expression error must improve on both stages, while fixed guards protect fidelity, reach, fold-change accuracy and discrimination; inspect the worst family as well as the mean. A failed confirmation must remain a failure. These guards are screening criteria, not a forecast of the official anchored score. A neutral-only success does not authorize claiming an improved submission.

## Reproduce from Windows Command Prompt

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe collect_post_v11_data.py
.\.competition_env\Scripts\python.exe test_post_v11_research.py
.\.competition_env\Scripts\python.exe audit_post_v11_transfer.py
.\.competition_env\Scripts\python.exe benchmark_post_v11_distributions.py
.\.competition_env\Scripts\python.exe audit_post_v11_composition.py
.\.competition_env\Scripts\python.exe document_post_v11_research.py
```

`plan_post_v11_experiments.py` was run once before future modeling. Its frozen JSON is tracked; rerunning it deliberately refuses to overwrite the plan. Neutral-test checkpoints verify plan and per-gene-file hashes before reuse. Data downloads verify cached bytes; changed upstream inputs require an explicit new snapshot rather than silently replacing this one. New code uses the existing competition environment and leaves its dependencies unchanged.

`test_post_v11_research.py` passes held-family leakage, masked-value abstention, no-source abstention, cosine missingness, count/depth-axis and valid-template checks. Documentation requires all 48 neutral records, exact-depth/count invariants and unchanged hashes for V11's numerical dependencies and 17 frozen generation snapshots. `research_manifest.json` records the complete input/output/code provenance. No submission was generated or uploaded by this work.
