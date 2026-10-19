# Research after V68: separate response information from cell generation

Review date: 19 October 2026. This is a literature and code review, not a new fitted model, experiment result or upload. V66 remains the highest published project score, **0.150882188**. The recommendations below are hypotheses; none establishes an improvement on the hidden competition panel.

## Recommendation

Run a small, crossed mean/decoder diagnostic first. If it shows that a more faithful count distribution recovers treatment effects at matched expression means, build a **control-anchored probability-distribution response predictor**. If accurate means rescue the existing generator while improved distributions do not, redirect effort toward a raw-count, uncertainty-aware response teacher and additional experimental contexts. Do not choose the next architecture from a paper's headline benchmark score.

We have tried many broad model classes. The remaining opportunity is more specific than another ensemble: learn the treatment-induced change in the probability of zero, positive expression and coordinated variation, without subsequently destroying those learned properties through exact count-margin repairs. This also requires better transfer validation. The full method is still untested here; some of its ingredients and motivating papers were already considered.

## What our evidence establishes

- V63/V66 improved expression prediction through assay-aware response transfer. V66's gain over V63 was only 0.000459393 overall, driven by reach; five other hosted components declined. Its small public count cohort predicted the wrong direction for five of six hosted components. More tuning on that cohort is weak evidence.
- V67's intervention-network correction reduced broad profile error by only 0.0018% and failed target confirmation. This fixed network model is not a promising release candidate.
- V68's intact-control population selection passed the complete fitted zero-response checks, but paired equal-study profile error increased **21.30%**. Two H1 targets with 157/137 actual DE calls had zero downstream calls in its balanced-real predictions. Baseline selection alone did not reproduce those induced responses. Only 4/10 fits reached the fixed convergence rule, so this is not a theorem about all population-selection methods.
- The same V68 null comparison found **545–5,328 significant genes** for the older synthetic zero-response generator across five public contexts, versus zero for the real-cell alternatives. This identifies a calibration problem on these cases, not the sole cause of the private score plateau or an estimated universal false-positive rate.
- V52's joint mean/rank head did not beat its controls. Split-view HCT116 mean/rank correlations were only 0.055/0.022. Learning noisy response summaries more aggressively can fit sampling error rather than biology.

These observations leave two unresolved problems: predicting the correct response and representing it as a realistic count population. Improving one can worsen the other. Public raw scores are not the leaderboard's anchor-scaled scores.

Local evidence: [V66](V66_LEADERBOARD_RESULTS.md), [V67](V67_RESULTS.md), [V68](V68_RESULTS.md), [joint mean/rank](V52_JOINT_RESULTS.md).

## What the literature adds

### 1. Test the output likelihood and sampling separately

[Bhattacharya et al., ExpressionVAE/design-space preprint](https://www.biorxiv.org/content/10.64898/2026.06.15.732063v1.full) compares discrete-latent priors and output heads. Its same-weight NB sampling/deterministic ablation changes distribution-sensitive and mean-sensitive metrics differently. This is evidence that decoder behavior can change conclusions without changing the trained predictor; it is not evidence that NB sampling improves all competition metrics.

The [author implementation](https://github.com/sanjukta7/sc-evae/blob/2b367c2acd085478705e2934edfdcf5c6121833b/src/sc-evae/models/output_heads.py) confirms that `NBHead.predict(..., sample=False)` returns its mean, while `sample=True` samples the same mean and learned gene dispersion. Its likelihood factors the mean into library scale and gene proportions. Dispersion is gene-specific and cell-independent, so it is not already a complete treatment-specific heterogeneity solution. Independent NB draws do not guarantee the specified library total exactly. The repository is MIT licensed.

Access limitation: the full bioRxiv page returned 403 and OpenReview requested browser verification. I inspected indexed primary-paper passages, the author project page and pinned decoder/VAE source. This review does not claim a full reading of an inaccessible PDF.

**Project implication:** isolate likelihood/sampling from learned mean, then score the final raw counts. Copying a new prior while retaining our existing emitter could miss the useful change.

### 2. Learn population distributions directly, but respect the task mismatch

[PopPert, September 2026 preprint](https://arxiv.org/html/2609.01357v1) predicts marginal-distribution residuals and low-rank Gaussian-copula dependence. Its main experiments use log-normalized zero-inflated Gaussian mixtures; raw-count ZINB is an alternative Norman experiment. Shared dependence loadings are fixed. Training uses pre-fitted target summaries and sampled pseudo-observations. Its HepG2 evaluation still includes 887 treated training conditions in HepG2; it does not demonstrate our controls-only unseen-context task. Main genetic outputs cover 2,000–5,000 selected genes, rather than all 18,533 competition genes.

The [author repository](https://github.com/whd1125/PopPert/tree/248fc3e1428de96335bef1d0271d30ca1e4d5565) was inspected at a fixed revision, including main and ZINB sampling code. Its README states that no software license has been selected. Use the scientific idea for an independent implementation rather than copying that implementation into the project.

**Not a newly discovered idea:** PopPert was cited in V17. V17 implemented a single shared random response gain, not the complete marginal/dependence predictor. The useful next distinction is direct count-distribution learning with an audited final decoder, not simply adding variance.

### 3. Fit sparse effects and their uncertainty from counts

[GPerturb, Nature Communications 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12215016/) uses hierarchical Bayesian sparse perturbation effects with Gaussian or count likelihoods. Its comparisons show that preprocessing/likelihood choices can substantially change inferred effects and their direction. The principal count version uses zero-inflated Poisson; the paper also considers Gamma-Poisson. This is not proof that NB is universally preferable. Its covariate/perturbation setup does not automatically solve anonymous unseen-context transfer.

The [MIT author implementation](https://github.com/hwxing3259/GPerturb/blob/05d5a431413ff9dc443c82b994b571ecc563426d/GPerturb/GPerturb_model.py) contains Gaussian, ZIP and ZINB variants. I inspected their model/loss interfaces. This review did not install or execute the remote package, and does not certify its training path.

**Project implication:** a second priority is a count-level source-effect teacher with uncertainty, compared against our existing inverse-variance pooling and guide-aware summaries. Replacing source summary labels with another neural network is not equivalent to fitting this likelihood.

### 4. Keep heterogeneous efficacy as a conditional hypothesis

[Song et al., Nature Cell Biology 2025](https://www.nature.com/articles/s41556-025-01626-9) studies continuous perturbation-response scores and nonlinear relationships with target suppression. Its PS implementation is available from the authors. Its genome-scale Jurkat experiment is high-MOI and has a targeted readout, making it unsuitable as a direct single-target full-transcriptome replacement dataset.

This was already researched in V17 and the post-V20 review. V17's response-gain test did not establish transferable benefit. Revisit only if raw source cells support reproducible treatment-specific mixtures after controlling for depth and batch. Hidden A/B/C efficacy cannot be estimated from their nonexistent observed treatment profiles.

### 5. Match the benchmark to the competition

[VCBench, April 2026 preprint](https://arxiv.org/html/2604.27646v1) separates unseen-context, unseen-perturbation and cross-dataset evaluation. It reports substantial sensitivity to splits, preprocessing and dataset aggregation. This supports testing the precise transfer task rather than borrowing performance claims from easier settings.

I also found a [public VCC competitor pipeline](https://github.com/rpinho/vcc2026), whose author reports 0.1528. Its transfer/spread-calibration approaches are close to our score range; it is not evidence of how the leading private models reach roughly 0.47. The competition permits undisclosed code/weights, so the leaders' architecture or data advantage cannot be inferred from rank alone.

## The next diagnostic: identify the bottleneck before training

This is a proposed design, not a frozen or executed protocol. Finalize sample eligibility, primary endpoints, seeds and thresholds before fitting or viewing new outcomes.

Use source populations with adequate controls and treatments. Split treated cells into a diagnostic fitting half and an independent scoring half, with disjoint NTC fit/reference cells. When guide/batch metadata permit it, split by those units, not just individual cells. Otherwise disclose shared-guide biology. Match population sizes between arms; do not manufacture 400 independent observations by duplicating a smaller treatment group.

Compare four arms for the same cases:

| Response-mean input | Count model | Question |
| --- | --- | --- |
| Frozen V63/V66 prediction, whole study excluded | Existing emitter | Reference |
| Diagnostic treated-half mean | Existing emitter | Does much more accurate response information rescue the old generator? |
| Frozen prediction | Distribution fitted from the diagnostic half, constrained to the supplied mean | Does better distribution shape help at the available predicted mean? |
| Diagnostic treated-half mean | Same fitted distribution family | Can this representation approximate the remaining treated half at all? |

The diagnostic arms intentionally use treated information and are **oracles**, never deployable A/B/C predictors or eligible release candidates. A fitted parametric distribution is not a true oracle: failures can reflect limited cells, family misspecification or fitting error. The two factors are not perfectly separable if shape constraints and means conflict. Report feasibility and realized means rather than claiming an orthogonal comparison when the final populations differ in both.

Score final normalized profiles, mean/bulk effects, detection, dispersion, rank effects, all available DE components, library-size distribution and null calibration. Report target panels large enough for discrimination; a two-target rank score is too coarse. Use matched shuffled shape descriptors and unchanged-control comparisons. Compare against the actual measured real-vs-real sampling floor. No mean of public raw metrics is called a hosted overall forecast.

If all four arms fail, the experiment has not established that another transfer network will help. If the fitted family cannot reconstruct an independent half of the same treated population, stop before cross-context training.

## Candidate if the decoder diagnostic succeeds

Independently implement a control-anchored integer-count probability transform. Learn source treatment residuals for zero/detection probability, positive-count shape and, only afterward, a small dependence component. Condition on query-control summaries and source evidence; predict parameter changes rather than cell-by-cell paired targets. Use raw counts with depth/batch adjustment. Compare a simple NB marginal against a hurdle extension before adding mixture complexity; many zeros alone do not establish excess zero inflation.

Use randomized ranks within discrete count masses and inverse count CDFs. When the response is zero and both CDFs are identical, the same uniforms must reconstruct the input counts exactly. Verify this identity numerically, including ties, clipping and large counts. A real baseline row provides dependence, while per-gene count changes can create perturbed states beyond V68's selection-only control-state hull. This preservation does not guarantee correct treatment-induced cross-gene dependence.

**Change the final constraints:** do not automatically force both exact V63 gene totals and exact original row libraries through the old dual-moment repair. Such repairs can change zeros and ranks. Conditional count generation can produce valid integer counts with fluctuating total libraries; validate their final normalized expression and competition size limits explicitly. If a mean constraint is necessary, fit it inside the probabilistic model and document any residual mismatch. Do not claim fixed-mean benefits when only the intended means match.

This differs from V15's five-quantile blend, V17's scalar response gain, V27's calibrated neural residuals and V34's flow-shape projection. Those tests did not establish that this specific direct count-CDF pipeline works; their failures also mean it needs negative controls, not an automatic upload.

## Validation and data priorities

1. After diagnostic reconstruction succeeds, train residual distribution parameters with complete study exclusion, including source response banks and preprocessing. HCT116/HEK293T form one Orion study; all related CD4 conditions belong together. Query NTCs are legitimate context input; query treated cells are not. Compare no-shape-change, source-pooled, shuffled-target and existing-generator controls.
2. Use target confirmation and study-level reporting, but acknowledge the existing biology has been reused extensively. More seeds are not more studies. A new donor within CD4 is additional donor evidence, not a new independent study; D3 remains unopened here.
3. Before a large new model, audit raw-data access for the previously identified KOLF/iPSC resources. The earlier KOLF audit found 273 panel targets, but only quantized mean summaries were used; those do not supply count distributions. The earlier full AnnData size was about 189 GB, so selective access and gene/guide/control coverage should be checked first. Sparse per-donor iPSC data can support fitting/uncertainty diagnostics without pretending to have 400 cells per target/donor. No new biological matrix was downloaded in this review.
4. If response means remain the dominant error in the diagnostic, prioritize a hierarchical raw-count source teacher and broader perturbation contexts over another generator. Use uncertainty and matched controls to test whether the teacher learns reproducible target effects beyond our existing pooling, rather than merely shrinking noisy labels.

The first milestone is an independently scored reconstruction/mean diagnostic, followed by a deployable whole-study test. This is a more informative use of compute than uploading another small blend with no confirmed public advantage.

## Reproduction and limits

[Root audit script](review_post_v68_literature.py) downloads ten author files at exact commits into an ignored inspection cache, records SHA-256 hashes and Python symbol inventories, binds twelve local evidence files, and independently recomputes the V68 equal-study profile ratio. Downloaded Python is parsed, never imported or executed. The compact receipt is [analysis/post_v68_research/source_audit.json](analysis/post_v68_research/source_audit.json).

```cmd
.\.competition_env\Scripts\python.exe review_post_v68_literature.py
.\.competition_env\Scripts\python.exe review_post_v68_literature.py --check
REM On a fresh clone, restore the ignored source cache first:
.\.competition_env\Scripts\python.exe review_post_v68_literature.py --restore-cache
```

The first command is for initial collection and refuses to overwrite its receipt. The second uses the saved cache without network access. An existing checkout without the ignored cache uses `--restore-cache` to reacquire and verify the exact pinned files without rewriting the receipt. Hash checks establish provenance/arithmetic, not model validity. This turn trained no new model, opened no reserve and generated no submission.
