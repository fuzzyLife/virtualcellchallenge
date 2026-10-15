# V26: learned populations with fixed V24 response means

Experiment started2026-10-15. V24 is the current published reference at0.14299785521745675. The user authorized the proposed count-decoder experiment. All prior numerical scripts, fits, plans, receipts and packages are preserved. No automatic upload or Git push.

## Why this is a distinct test

V17 sampled scalar strengths; V21 modified limited pretrained residual programs; V22 fitted bounded linear state-response interactions inside the same expected/Poisson/empirical generator. V23 repaired its emitted means. V26 replaces the trainable distribution mechanism with a neural logistic-normal multinomial decoder trained on actual counts. It starts by testing real-control calibration, which returning V24 bitwise cannot establish. V22's factorial oracle already supplied the reason to prioritize distribution error; it is not repeated or described as a new result.

## Trainable model and data

The posterior uses4096 count-residual features, eight stable-state factors and16 response descriptors, with a128-unit hidden layer and16 stochastic latent dimensions. The decoder uses latent/state/response inputs,64/128 hidden units and18533 output logits bounded+/-4. Counts have a multinomial likelihood conditional on observed aligned-axis library sizes. This is a learned correlated logistic-normal mixture, not a claim that independent Poisson counts reproduce biology.

State features are512 training-only, shared nonpanel genes with relatively stable source responses; an eight-dimensional control PCA supplies state coordinates. The response descriptor reuses the first16 directions of V25's frozen training-only intervention basis. It does **not** reuse V25's rejected covariance transfer operator. Stable-state measurements in treated cells remain imperfect basal proxies; source paired cells or individual counterfactuals are not observed.

Train three matched networks: independent state input zero; conditional state input; state-shuffled training with state vectors permuted within source batch and fit/validation partition. Source-batch training shuffling disrupts state/count dependence; permuting generated cells alone would not test it. All networks use the same samples, architecture, initialization,12 fixed epochs, AdamW learning rate0.001 and KL coefficient0.001. Choose epoch by held-cell validation ELBO only. This split does not make an independent biological study.

Source contexts are H1, HCT116, HEK293T, resting CD4D2 and stimulated48h CD4D2. Exclude the complete evaluated family from raw-count training, feature selection, response basis and V24 response estimation. Excluding CD4 excludes both conditions and all related mean specialists. Donor aliases across summary/raw releases are not proven identical.

Sample at most64 controls per source sequencing batch, and96 treated cells per source/target. Require adequate control support and at least16 treated cells. Split fitting/validation before estimating group means. Fit-only means receive0.001 control-profile shrinkage and1e-8 measured-gene pseudomass to avoid impossible likelihoods for rare validation-only detections. This is likelihood regularization, not the rank-correlation add1 operation previously avoided. Original row IDs, labels, batches, fit flags, aligned libraries and complete-source library metadata are saved. The aligned-axis sums are appropriate for this output-count likelihood; they do not replace complete-source normalization in earlier models.

The18533-gene output head learns residuals only on genes measured in every training source, excluding the300 target-panel genes from learned state/response residuals. Unsupported coordinates retain empirical query variation; the direct target retains the reference mean convention. Coverage is reported. This is a full output-axis count generator with an explicit learned-coordinate mask, not a claim of learned effects for unmeasured genes.

Scientific runtime: unchanged Python3.12.10, NumPy2.5.4, SciPy1.18.1, pandas3.0.6, PyTorch2.10.0+cu126, RTX2070SUPER8GiB. Two BLAS/CPU threads, deterministic Torch algorithms, fixed seeds, no TF32. Saved weights have no pickle objects. Emission uses fixed NumPy latent draws and CPU inference so GPU scheduling cannot change packaged counts. Official packaging would use the separate competition environment.

## Mean and library contract

For every paired case, generate the family-excluded V24 comparator first. All new count arms preserve its exact cell libraries and gene totals. They also limit mean per-cell CPM probability L1 drift to0.0005. No real treatment means enter inference or repair.

Generate conditional gene probabilities, fit mean/pooled proportions by damped exponential tilting, draw multinomial counts, repair exact gene totals within16 library-depth strata, and perform bounded2x2 integer swaps to improve the per-cell mean constraint. Final swaps preserve global gene totals and cell libraries; they can alter depth-stratum margins. The CPM-swap budget is5% of library molecules and is logged. These operations may change learned shape, so all treatment comparisons evaluate the final integer counts. The empirical comparator receives the same sampling/repair machinery.

Before freezing, synthetic400x18533 tests exposed mean drift0.00511 with pooled totals alone and0.00233 with depth-stratified totals. The final bounded swap solver reduced the latter to0.000000975, preserving both exact margins. No biological treatment score was opened during this numerical design. Four focused tests verify gradients, CPU decoder agreement, feasible moments, deterministic integer margins, zero support and swap invariants.

## Prospective calibration and release gates

Six arms: V24, empirical-query distribution, independent neural, conditional neural, target-descriptor-shuffled neural, and state-shuffled-training neural. Only conditional is an eligible finalist; diagnostic arms cannot be promoted after inspecting results.

For each context, combine the old left/right control pools, then freeze a new three-way split:400 actual null cells,400 reference cells, and remaining inference controls containing400 distinct template cells. This corrected the unsupported attempt to draw400+400 from H1's768-cell old left half; that preparation aborted before writing a plan. No cells are duplicated, and treatment splits remain unchanged.

Train the three H1-excluded models, evaluate its null gate, and stop if it fails. Otherwise continue HCT116, HEK293T and both CD4 states. Conditional must satisfy the mean contract and all of these per-context real-null envelopes:

- DE calls <=max(20, real-null calls+10).
- Detection MAE <=max(0.03,2*real-null detection MAE+0.005).
- Log-variance MAE <=max(0.5,1.5*real-null log-variance MAE+0.1).

These are prospective engineering gates, not calibrated universal confidence limits. One case per context is an initial falsification screen; passing it alone cannot establish distributional generalization. Failure halts before treatment outcomes and preserves the existing reserve. This new protocol does not relax any old rejection policy.

Only after all null gates pass, execute the44 reused public treatment cases/six arms. CD4RestD1 has200 real treatment cells; the other36 cases have400. HCT116 NMAE remains missing. The baseline is the V24 algorithm under the same source exclusions, not the published ABC package. Keep old available raw guard thresholds, >=3 improved primary metrics including DE, worst-family fidelity/reach, the mean contract and improved centered shape distance. Conditional must also beat empirical, independent, target-shuffled and state-shuffled training on at least two DE metrics each.

Only one locked passing finalist may open the new raw CD4D3Rest reserve, excluding all CD4 sources. That reserve's file is still partial from V25 and must first be completed, author-MD5 verified, extracted and independently row-audited. No data or validation success is claimed in advance. A release additionally requires whole-package count/axis/provenance and exact replay checks. Failure retains V24 without creating a duplicate upload.

## Reproduce

Run in Windows Command Prompt using the pinned scientific environment. Frozen plans/decision files refuse overwrite; identical verified training caches can resume.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
set "VCC_SCI_PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
"%VCC_SCI_PY%" test_v26_decoder.py
"%VCC_SCI_PY%" prepare_v26_experiment.py
"%VCC_SCI_PY%" run_v26_with_recovery.py
"%VCC_SCI_PY%" audit_v26_experiment.py
```

Root scripts and small plan/model/score/hash evidence are versioned; raw sparse training subsets, model arrays and population caches are ignored. V26_RESULTS.md and V26_SUBMISSION.md record actual execution and the artifact decision after the gates finish. Historical best statements in old documents remain snapshots of their publication time.

## Execution reporting recovery and descriptive diagnostics

The first complete null assessment encountered an operational type error while hashing its filename: the existing sha256 helper expects a Path, while the frozen assessor passed a string. Counts and checks had completed, but no decision file was written. assess_v26_recovery.py runs the identical hash-pinned assessor with Path coercion only; run_v26_with_recovery.py dispatches that adapter from the identical frozen pipeline for from-scratch reproduction. No model, output counts, samples, thresholds or decision logic changed. analysis/assessment_recovery.json preserves input table hashes, original assessor hash and adapter hash. The original pipeline then resumed verified caches and saved its planned stop.

Source latent diagnostics use1276 held-validation control cells: conditional aggregated posterior variances range0.919..1.062, median0.992; its mean-vector norm is0.102. These diagonal moments do not show gross source-prior variance collapse, but do not establish full joint matching or unseen-context calibration. Correlations and a diagonal Gaussian moment proxy are descriptive, not exact mixture KLs or causal estimates. Root diagnose_v26_latents.py records the hashes and measurements without changing predictions or gates.

Prerequisites include the verified raw extractions, source_cache, V24 teacher and V25 held-family intervention bases referenced by the plan. Reproduce those from the [V24](V24_STRATEGY.md) / [V25](V25_STRATEGY.md) strategy documents before V26 on a fresh checkout. The scientific environment also requires Torch2.10.0+cu126 alongside requirements-v23-science.txt; current actual versions are archived in analysis/runtime.json.

After calibration rejection, diagnose_v26_null_genes.py exactly replays the saved conditional null and classifies errors by learned/fallback mask and observed source support. It is descriptive post-rejection analysis of spent development biology, not new independent validation or a model amendment. Its category table, per-gene calls and hashes are retained; no revised mask or sampler was deployed in V26.
