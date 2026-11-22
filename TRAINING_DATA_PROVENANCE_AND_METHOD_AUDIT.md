# Training-data provenance, model pipeline, and project findings

**Project:** Virtual Cell Challenge 2026
**Updated:** 22 November 2026
**Scope:** Development through the final V91 artifact; this is a provenance and methods record, not a claim of final-test accuracy.

This audit explains which data actually trained or informed our perturbation predictions, how we obtained and transformed it, which code used it, what we tried but did not put into the final predictor, and what we learned by comparing our work with the public CellularAI release. “Training data” here means measured perturbation-response evidence used for transfer. Challenge non-targeting controls are baseline populations, not treated training examples.

## 1. What the competition asked us to predict

For each of three held-out final contexts D, E, and F, the challenge supplied 18,400 non-targeting control cells and a list of 300 CRISPRi target genes. We had to predict 400 integer-count cells for each target and context on the 18,533-gene challenge axis: 300 × 3 × 400 = 360,000 cells in 900 groups. The validation phase used contexts A/B/C; final contexts must be generated from their own controls and must not be created by relabeling A/B/C predictions.

The final controls tell us each context’s baseline expression, cell-to-cell variation, and library sizes. They contain no perturbation outcomes. The perturbation-response estimates therefore came from external public studies, transferred and calibrated against those baselines. We did not have final D/E/F treatment labels.

## 2. Perturbation-response data actually used

The frozen final response cache is audited in **experiments/final_v91/analysis/source_audit.json**, **source_provenance.json**, **target_coverage.csv**, and **source_cache.npz**. It has 1,487 target-by-source response records: 638 carried forward bitwise from the earlier cache and 849 newly extracted. Its five source families are K562 (260 rows), CD4 (625), H1 (28), HCT116 (282), and HEK293T (292). These are summary rows, not 1,487 independent experiments or cell lines. CD4 includes three correlated conditions. In total, 298 of 300 final targets have at least one eligible response. ANG and DPY30 have none; their downstream response stays at the control-based baseline, subject to the fixed on-target assumption.

| Source | Where it came from and how we collected it | Transformation and use |
|---|---|---|
| K562 | The project’s public K562 GWP screen, stored locally as **reference/K562_gwps_raw_bulk_01.h5ad**. Its checksum is checked against the earlier source receipt. Acquisition/preparation began with **download_reference.py** and **build_submission.py**; final-panel extension is explicit in **build_final_sources.py**. | The final builder parses target and non-targeting labels, uses eligible guide rows and source controls, requires at least 30 cells, aggregates counts weighted by cell number, maps source symbols to the challenge gene axis, and computes log2 CPM response summaries, standard errors and adjusted p-values. The target coordinate is masked so a guide’s direct knockdown is not learned as a downstream response. |
| CD4 | Public GWCD4i differential-expression summary objects, including a merged summary and two donor-condition groups, in **GWCD4i.DE_stats.h5ad/.h5mu**. **extract_cd4_reference.py** pins the S3 object identities and ETags and reads needed byte ranges. **collect_final_cd4.py** extracts the 171 final targets absent from the old panel and verifies old/new row identities and axes. | **build_cd4_model.py** applies eligibility and quality filters to log fold change, adjusted p-value and standard-error fields. The final source builder preserves reused validation-target rows and adds eligible new-target rows. The three conditions are correlated evidence, not three independent cell lines. |
| H1 | Arc Institute’s public 2025 Virtual Cell Challenge H1 train, validation, and test H5AD files, accessed by pinned generation/ETag and byte-range reads. Identities are recorded in **experiments/final_v91/analysis/H1_new_source.json** and **extract_v5_h1.py**. **collect_final_h1.py** selects newly needed targets; shared non-targeting controls are retained once across splits. | Selected raw counts are mapped by gene symbol to the challenge axis and summarized with the source-response routine. The final extension comprised 8,463 treated cells for 12 new targets. A target/source summary must meet the minimum cell requirement to enter the response cache. |
| HCT116 and HEK293T | Public X-Atlas-Orion release in Hugging Face dataset **slaf-project/X-Atlas-Orion**, pinned at revision **598aa5442a4ee9a8f1a383a843ab5ce142904d04**. The original stratified selector is **prepare_v7_sources.py**; **collect_final_orion.py** extends the same frozen source to new final targets, using the reviewed extractor and offset repair. | The selector uses full source metadata, filters to passing guide/cell records, samples up to 400 cells per target across disjoint batch partitions (200 per partition), and samples up to 32 non-targeting controls per batch. For the final extension, the pinned release included 30,795 HCT116 cells across 109 batches and 42,155 HEK293T cells across 223 batches, including controls. Extraction verifies cell IDs, gene maps, offsets, full library totals and hashes before pseudobulk summary. Missing measured genes remain masked, not set to biological zero. |

The final cache has uneven coverage: 19 targets have five source families, 180 have four, 82 have three, 14 have two, three have one, and two have none. The source audit separately records the 129 final targets shared with validation and 171 new targets. The 638 reused records are not new replication. Evidence from a single source is not equivalent to evidence from several independent lineages.

### Collection and transformation code map

The end-to-end final acquisition and model path is recorded in **PROJECT_PROGRESS.md**, **FINAL_ROUND.md**, and the receipts under **experiments/final_v91/analysis/**. Main code:

- **download_reference.py**, **build_submission.py**: initial K562 source acquisition, summary construction, gene/context normalization, common file-contract utilities, and early baselines.
- **extract_cd4_reference.py**, **build_cd4_model.py**, **collect_final_cd4.py**: pin and extract CD4 DE evidence, apply source filters, append only new final targets.
- **prepare_reference_controls.py**, **extract_v5_h1.py**, **collect_final_h1.py**: use pinned H1 byte ranges, select cells from the public 2025 splits, verify shared controls and summarize responses.
- **prepare_v7_sources.py**, **collect_final_orion.py**: pin and sample X-Atlas-Orion HCT116/HEK293T sources; final extraction preserves source IDs and offsets.
- **fit_v7_response.py**, **fit_lineage_submission.py**, **build_final_sources.py**, **audit_final_sources.py**: construct the normalized target/source response cache, retain standard errors, significance and knockdown fields, preserve missingness, exclude each perturbation’s own coordinate, and audit coverage and hashes.
- **inspect_final_inputs.py**, **diagnose_final_contexts.py**: validate final controls, gene axis and target manifest; calculate descriptive control-profile similarities. Similarity is not a cell-line identity or response-accuracy result.
- **final_response_methods.py**, **fit_final_submission.py**, **audit_final_model.py**: fit selected response corrections and priors, freeze inputs/code, and check numerical bounds and identity guards.
- **build_final_submission.py**, **summarize_final_generation.py**, **package_final_submission.py**, **validate_final_submission.py**, **run_final_pipeline.py**: synthesize counts from D/E/F controls, stream and package the 360,000 cells, run official-format and independent full-payload validation. The orchestrator does not upload or push.

## 3. How the final V91 predictions were formed

The final D/E/F recipe is a refit of the best published validation approach, V90, on the actual final targets and final controls. It is not a model trained on hidden final outcomes.

1. The source cache retains per-gene response estimates, standard errors, adjusted p-values, knockdown information, measured-coordinate masks and source-family labels.
2. **final_response_methods.py** and **fit_final_submission.py** apply the empirically retained V13/V16 source pooling, shrinkage and missingness handling.
3. They apply the V24 **cdf50_de** donor-supported CD4 correction; the V66 **assay75** strength correction under the existing 80% V24 identity-margin guard; and the V90 **interval25** normal-moment interval projection at quarter strength, capped at 0.125 LFC, under the existing 80% V66 identity guard.
4. The chosen response profile is sent to the unchanged **residual10** count emitter using actual final D/E/F control cells, their control distributions, and their library totals. Random seed is 20261011. The direct on-target factor is 0.1 and is applied once. The model uses no final treatment measurements.
5. The official CLI 0.3.0 pack and independent checks verified final context labels, 300 targets per context, 400 cells per group, 18,533 genes, all 360,000 cells, library totals, cell identities, sparse payload and regenerated groups.

The exact artifact and checksums are in **FINAL_V91_RESULTS.md** and **experiments/final_v91/analysis/submission_ready.json**. Validation establishes a technically valid and reproducible file; it does not establish biological accuracy or a better final rank.

### Results from the control-cell and neural-network analyses

These analyses helped characterize the query contexts, but they are not supervised perturbation training labels.

- **All 55,200 validation controls were analyzed** (18,400 each from A/B/C). In the full measured-gene space, each context formed its own pure cluster: an unsupervised K=3 rank-correlation clustering assigned every A, B, and C cell to separate clusters (purity and ARI 1.0). This confirms strong context structure in the measured matrix; it does not establish biological identities or response-transfer accuracy. See [full-matrix clustering results](FULL_RANK_CLUSTERING_RESULTS.md).
- **External control similarity was broad and depth-sensitive.** Per-cell full-gene Pearson medians were A: CD4 0.482, H1 0.381; B: H1 0.464, CD4 0.400; C: H1 0.477, CD4 0.415. The closest external-control assignment was essentially all CD4 for A and H1 for B/C, but these are profile similarities, not cell-line labels. The sampled HCT116/HEK293T controls in this diagnostic came only from batch 1. See [reference correlations](CELL_REFERENCE_CORRELATION_RESULTS.md).
- **Restricting to exactly the 300 perturbation target genes changed the correlation view.** Pairwise-valid-value Spearman ranks used the original counts, average ranks for ties, and no pseudocount, transformation or imputation. Median A correlations were CD4 0.513 and H1 0.391; B was H1 0.463 and CD4 0.414; C was CD4 0.416 and H1 0.402. K562 covered only 231 target genes; other references covered all 300. Pooled A/B/C correlations were 0.628 for A–B, 0.525 for A–C, and 0.656 for B–C. See [target-gene rank results](TARGET_RANK_CORRELATION_RESULTS.md).
- **Cell-to-cell spread was substantial.** For the 300-gene profiles, individual-cell correlations varied widely within each context. Own-context leave-one-out profile correlations were around 0.67 on average; sampled same-context cell-pair Spearman medians were approximately 0.48–0.51. Correlations to reference profiles were also strongly associated with total counts and detected target genes. These diagnostics cannot separate biological variability from capture depth or dropout. See [target-cell dispersion results](TARGET_CELL_DISPERSION_RESULTS.md).
- **Neural nets reconstructed control expression, not CRISPRi effects.** We trained 24 context-specific candidates (four architectures, normalized/raw objectives, three contexts), selecting on seven held-out guide groups. Best normalized networks beat the mean baseline by about 8.0–10.9%, while beating the ridge model by only 0.85–1.30%. Raw-count networks improved 26.0–31.9% over a depth baseline, but that objective was not directly comparable to normalized-profile loss. This showed the 300-gene input carried some predictive signal, but did not show causal perturbation learning and did not replace the response-transfer model. See [neural-network results](GENE300_NEURAL_RESULTS.md).

For K562 specifically, the final source builder keeps control rows marked non-targeting, core-control, and at least 30 cells. A treated guide row is eligible when it has at least 30 cells, finite fold-expression, and fold-expression at most 0.5. It aggregates source counts weighted by cell number, converts treatment and control pseudobulks to CPM, computes log2((treated CPM + 1)/(control CPM + 1)) clipped to [-3, 3], estimates standard errors from the two pseudobulk counts, obtains normal-approximation p-values and BH-adjusts them over genes with control CPM above 5. The target’s own coordinate is excluded from downstream validity. These details live in **build_final_sources.py** and the historical K562 model receipt.

## 4. Other data we gathered, and whether they entered the final model

A recurring source of confusion was treating every downloaded dataset or model as “training data.” These experiments had different roles:

| Data or project | Collection and analysis | Did it feed final V91 response predictions? |
|---|---|---|
| Challenge controls and specifications | Downloaded challenge manifests, non-targeting profiles, target/cell counts, and the 18,533-gene axis. We analyzed cell-type correlations, target-gene correlations, per-cell dispersion, clustering, and context mapping. Scripts include **analyze_cell_reference_correlations.py**, **analyze_target_rank_correlations.py**, **analyze_target_cell_dispersion.py**. | Final D/E/F controls fed the baseline and count emitter. Validation controls informed diagnostics and selection, not final treatment truth. No validation or final perturbation truth was available to train on. |
| RPE1 public study | **collect_v13_rpe1.py** retrieved Figshare article 20029387, file 35775581. | Exploratory transfer tests only; not one of the five frozen final response families. |
| PIE | The public PIE release was cloned/reviewed; the cross-dataset checkpoint and saved training-response evidence were audited. **infer_v18_pie.py** and V56/V58 experiment code tested all three outputs, separated 154 seen and 146 unseen targets, compared against shuffled/evidence-only controls, and tested response prediction before count generation. Saved evidence contained 64,260 context-perturbation records from Tahoe, Jiang, VCC25 and Orion. | No. PIE was not added to the final response cache. Its means sometimes offered small response-level gains, while joint outputs and count decoders did not robustly beat V24; tests did not establish a winning replacement. See **V56_PIE_RELEASE_REVIEW.md** and **V58_PIE_CROSS_RESULTS.md**. |
| KOLF | **collect_v17_public_tables.py** inspected public metadata. The roughly 189 GB raw dataset was not downloaded and no response model was trained on it. | No. |
| LINCS fingerprints and network features | Public drug/perturbation fingerprints and interaction features were evaluated as potential transferable features. | Exploratory only; no direct final response labels. |
| Arc essential screens | **V83_COVERAGE_AUDIT.md** reports full target-vocabulary checks against the challenge panel. | No. Relevant essential-screen tables had zero overlap with challenge targets, so they supplied no response labels. |
| BOLERO | Code and paper were reviewed and a compact adaptation was explored. | No full BOLERO checkpoint was downloaded or used as the final predictor. See **BOLERO_REVIEW_V87.md**. |
| CellularAI / Billion Cell Atlas | Public repository, method guide, example config, release article and local clone were reviewed. | Code and data were comparison evidence, not our model inputs. We did not have Illumina’s BCA response bundle. See **CELLULARAI_REVIEW.md**. |
| Claude gap-analysis script | The user-provided **vcc_gap_analysis.py** was read. | No. It was a synthetic toy analysis, not a fit on our real cache or a source of model parameters. |

## 5. What the CellularAI comparison clarified

The [official CellularAI method guide](https://github.com/Illumina/CellularAI/blob/main/docs/method.md), [public example configuration](https://raw.githubusercontent.com/Illumina/CellularAI/main/configs/example.yaml), [Illumina article](https://www.illumina.com/science/genomics-research/articles/introducing-the-cellularai-algorithm--predicting-single-cell-res.html), and [Arc Institute challenge wrap-up](https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up) show that the first-place approach was a complete evidence-to-cells pipeline, not merely a larger neural network.

| CellularAI stage | What its public release describes | Our difference |
|---|---|---|
| Response evidence pooling | Per-gene random-effects DerSimonian–Laird meta-analysis, inverse-variance weighting, knockdown efficiency, query-control similarity, noise reliability, and paired/auxiliary-record factors. A separate weighted Stouffer/Liptak combination uses signed evidence. | We pooled a much smaller and more heterogeneous five-family bank with empirically tested family weights, shrinkage, missingness handling and guarded corrections. We did not reproduce the complete gene-wise random-effects and signed-evidence pipeline at broad lineage scale. |
| Context matching and adaptation | Control-profile similarity (the released example has sim 0.25 and sim_hvg 2000), plus a rank-160 response-residual basis and gene-wise ridge regressions from source controls to response coordinates. | We measured context correlations and tried lineage/donor matching and response priors. We did not have enough independent response-profiled lines to fit the same broad low-rank context correction. Correlated CD4 conditions did not make independent contexts. |
| Knockdown, DE plan and power | Uses measured knockdown information and an explicit power model; the released example lists kd_ref 0.6, power_k 111 and power_depth 20000. Plans which genes should respond using significance, LFC, query CPM, response breadth and detectable power. | We calibrated response strength and on-target handling, but our selected final method did not contain an equally explicit power-aware gene-call planner. Early conservative signal thresholds missed many response genes. |
| Cell synthesis | Starts from actual controls and targets per-gene mean/DE behavior using stratified CPM quantiles, planned Wilcoxon z scores, tie-corrected bisection, a Gaussian copula and repair/gating. The example z_plan range is 5.5–9. | We also generated cells from real controls and preserved library totals, but our residual10 emitter did not solve the same per-gene Wilcoxon-z targets. Matching total counts or moments is insufficient for rank/significance metrics. |
| Evidence scale | Illumina reported access to the Billion Cell Atlas and far broader cell-line coverage; the public release omits the private BCA-derived response bundle. | Our audited bank has five source families, 1,487 response records and eligible support for 298/300 targets. We could not recreate their private Atlas evidence from public code alone. |

These settings are examples from the public release, not proof of the private winning run’s exact configuration. The method guide also says raw preprocessing and generation of the summary bundle are outside the public repository. Cloning their code could not recover private measurements. We did not miss the general idea of transferring measured perturbation responses: V4 through V90 already did that. The gap was evidence breadth and independent context coverage, plus more direct planning and calibration for the metrics.

### What the screenshot leaderboard taught us

A “mean predictor” predicts the expected or pseudobulk expression change for each gene after a perturbation; it does not predict which individual cell will express which counts. Such means can be strong because many perturbations have reproducible average responses, and the challenge’s expression and DE metrics reward getting these average shifts and directions right. A single mean profile still cannot satisfy the required 400-cell output or reproduce cell-to-cell dispersion; it must be paired with a calibrated population generator.

The leaderboard image shared in this discussion showed leading entries with strong scores across several scaled metrics, not a single standout score on one metric. The challenge defines zero as the context-mean baseline and one as replicate-level performance. It averages six complementary scaled metrics: perturbation discrimination, expression MSE, DE log-fold-change NMAE, DE direction fidelity, DE direction reach, and DE significance overlap. Therefore a good model needs perturbation discrimination and expression accuracy, while also recovering DE magnitude, direction, reach and significance. Optimizing one proxy or synthesizing more variable cells can damage another metric. Rank also moves as other teams submit, so hosted Overall is more interpretable than rank alone.

## 6. Full published validation record and what improved

These are all 22 published validation scores supplied in this project history. Rank varied over time; Overall is the comparison recorded here. V91 is not a validation submission and has no known final-test result in the project records.

| Version | Overall | What the experiment showed |
|---|---:|---|
| V1 | -0.276290 | Initial public-response transfer was far below the context-mean baseline. |
| V2 | -0.223325 | Donor-calibrated hybrid improved, but remained weak. |
| V4 | -0.050674 | Lineage-guided response coverage made a large early gain. |
| V5 | -0.113637 | Atlas-inspired calibration regressed. |
| V6 | -0.056790 | Fidelity guard did not beat V4. |
| V7 | -0.051387 | Resampling was close to V4; effect release was essentially a fallback. |
| V9 | -0.070447 | Activation-preserving calibration regressed. |
| V10 | 0.116754 | Pooled intervention population transport was the first major improvement. |
| V11 | 0.137004 | Consensus-supported coverage and dispersion improved further. |
| V13 | 0.137883 | Assay-aware response transport gave a small gain. |
| V16 | 0.141614 | Residual/dispersion correction gave another small gain. |
| V20 | 0.139012 | Stacking regressed slightly versus V16. |
| V24 | 0.142998 | Context-matched donor response became the leading project score. |
| V36 | 0.141845 | Retrieval-contrast did not beat V24. |
| V50 | 0.141523 | Empirical response prior did not beat V24. |
| V55 | 0.138767 | Fixed V24/V36 half-blend regressed. |
| V61 | 0.116763 | PIE/state-response adaptation regressed. |
| V63 | 0.150423 | Identity-guarded assay response improved strongly. |
| V66 | 0.150882 | Assay-strength calibration gave a small additional gain. |
| V69 | 0.095263 | Probability-coupled count mixture regressed sharply; count realism alone did not solve the score. |
| V82 | 0.150882 | Predictions were identical to V66 and score matched; this was a reproduction/file-identity check, not an independent method win. |
| V90 | 0.152144 | Conservative interval projection reached the best published validation score in the recorded history, only 0.001262 above V66. |

The progression shows that pooled response transport and guarded, modest corrections helped on this validation panel. Several elaborate models, ensembles, PIE blends, and count-distribution changes failed to improve. This does not establish the same ordering on the hidden final panel. The V90 gain was small; repeated reuse of public cases limits the independence of our model-selection evidence.

## 7. Audit of the Claude script shared in the discussion

The supplied **vcc_gap_analysis.py** fabricates a 50-gene dataset, four source response matrices and five perturbations, then declares HCT116 to be “truth” and the strongest/closest source. Its reported LFC shrinkage, DE call counts and percentage gains are outputs of that synthetic setup, not measurements from our real data or historical submissions. Its ideas about source similarity, evidence quality, statistical power and metric-aware synthesis are reasonable hypotheses and overlap with CellularAI’s public method, but the claimed “31%” gain and response magnitudes must not be attributed to our model.

There is also a sign-dependent bisection issue: the script updates lower/upper bounds as though its rank statistic always increases with a gene shift. For negative target shifts, the statistic moves in the opposite direction, so its update can move the candidate the wrong way. A production Wilcoxon-targeting solver needs separate positive/negative tests, tie handling, zero-effect checks and comparison to an independent rank-statistic implementation.

## 8. Reproduction, provenance and limits

- Start at **PROJECT_PROGRESS.md** for the complete version history and source map; **FINAL_ROUND.md** and **FINAL_V91_RESULTS.md** record final-phase preparation, artifact checksum and commands. **CELLULARAI_REVIEW.md** gives the detailed code-level winner comparison. PIE, BOLERO, essential-screen and earlier experiment documents named above preserve individual tests and negative results.
- Final source receipts and compact audits are under **experiments/final_v91/analysis/**. They pin hashes and record source IDs, selected rows, gene mapping, coverage, missingness, response construction, final checks and code hashes. Large source matrices and challenge payloads are local/ignored; a Git clone alone is not a copy of the data.
- Reproduce with the project environment and documented pinned source identities. The final orchestration sequence is in **FINAL_ROUND.md**; it does not upload predictions.
- Do not turn absent response measurements into evidence for a biological zero. ANG and DPY30 remain unsupported and masked in the response source cache.
- The project record last reports a user-reported V91 upload in progress on 23 October 2026. It contains no verified server receipt, entry ID or final score. This document does not claim that the file was accepted or scored. The challenge deadline has passed; check the account/CLI receipt before citing a final-phase outcome.
- The official [evaluation page](https://virtualcellchallenge.org/evaluation) and [FAQ](https://virtualcellchallenge.org/faq) describe payload and scoring rules. Follow the organizer’s current records if they differ from this historical project log.
