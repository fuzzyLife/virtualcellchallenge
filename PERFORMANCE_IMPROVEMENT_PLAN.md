# Performance improvement plan

Research checked: 2026-10-07. This records the original proposed sequence preceding V7; the completed experiments and published result are recorded below.

## What our evidence says

V4 remains our best published candidate: overall -0.050674, versus V5 -0.113637 and V6 -0.056790 on the recorded validation panel. V4's scaled fidelity is -0.683855 and significant-set overlap is -0.029878. Its scaled expression score of zero is a clamped score, not zero error. These aggregate results cannot identify which context or target is responsible. Fidelity yield is not simply sign precision.

The control-expression neural networks learn observational relationships. A good prediction of other genes from the 300 measured genes does not establish their response to CRISPRi. The later embedding pilot does train on interventions: 889 source contrasts across 284 targets. However, it uses a single ridge penalty of 10, fixed response-label shrinkage and a fixed 25% blend. In its H1 public-source generated-cell test, Tahoe/rank reduces fidelity from 0.145431 to 0.057594, overlap from 0.014226 to 0.008027 and worsens the corrected expression ratio from 1.276690 to 1.283540. These are raw public development metrics, not leaderboard scores. Very weak predicted responses make underfitting a plausible explanation to test; this does not prove that lowering regularization will improve transfer.

Evidence and reproduction: [embedding pilot](EMBEDDING_RESPONSE_PILOT.md), [upstream review](UPSTREAM_REPOSITORY_REVIEW.md), [V6 results](V6_LEADERBOARD_RESULTS.md).

## Public approaches and material gaps

### Broader measured CRISPRi coverage

[AtlasShift](https://github.com/kaipengm2/Virtual-Cell-Challenge-2026) pools full prepared K562, HCT116, HEK293T, H1 and CD4 statistics, then reanchors responses on target-context controls. Its historical September 10 score is author-reported, not a current ranking or independently replicated result. Its preparation requires substantial downloads/scans.

Our K562 pseudobulk aligns only 7,681 of the 18,533 required output genes, despite covering 260 target genes. At the time this plan was proposed, HCT116/HEK293T were pilot batches. V7 subsequently extracted panel targets across their releases and retained 282/291 eligible contrasts, with stratified selected cells rather than every raw cell. Expanded estimators failed local guards and the published V7 reverted to exactly V4 effects. H1 covers 25 targets. Availability masks remain necessary; control-marker similarity does not establish an A/B/C intervention identity.

### Response pooling and uncertainty calibration

[Sidechain's model log](https://saberhq.com/sidechain/models/) describes combining effects before adaptive shrinkage, source-confidence floors, measured knockdown calibration, and controlled count-emission experiments. It reports SER-14 at 0.1626 on October 3; its Tahoe-neighbor variant changed 0.1131 to 0.1132. These are participant-reported results and do not establish causality or guarantee gains here.

Our CD4 estimator already shrinks effects using standard errors, and K562 uses a noise threshold. The missing experiment is learning the shrinkage scale for a target's combined response and accounting for disagreement among sources, rather than applying only fixed source rules. Count-generation calibration is also not new to us: V5 tried separate mean-CPM and pooled-profile matching, alongside other changes, and scored worse. It therefore needs an isolated test, not automatic reuse.

### Perturbation-conditioned neural models

[Arc STATE](https://github.com/ArcInstitute/state) distinguishes a cell embedding model from a transition model trained to predict interventions. Its repository provides Replogle-Nadig genetic-perturbation training and whole-context zero-shot split examples. This is a relevant neural baseline, not proof of a superior current private leaderboard entry. Any comparison needs genetic-perturbation training, complete output-axis handling and conversion to the competition's raw count contract.

The [Stage 1 Lab](https://github.com/danielchen26/vcc2026-stage1-lab) explores rank-test detectability and count construction, but its README marks the biological response estimator unfinished. Treat its scorer claims as hypotheses to verify against all six official metrics, not as a demonstrated winning solution.

## Controlled sequence

1. **Expand source coverage first.** Cover the complete HCT116 and HEK293T release with stratified panel-target and matched control cells, and inspect the full K562 raw axis before deciding whether it expands coverage. The V7 probe subsequently confirmed that the full GWPS raw file has the same 8,248 measured genes (7,681 challenge overlaps); it cannot expand the output axis. Record measured-gene masks, targets, guides, batches, donors, knockdown depths and uncertainty. Deduplicate controls across public splits. Freeze versions/hashes. Cache summaries so fitting does not require repeated full atlas reads.
2. **Create stronger transfer baselines.** Compare the original V4 estimator with full-source pooling, confidence floors and target-adaptive shrinkage. Test normalization by measured knockdown strength with conservative caps; this assumes transferable dose dependence and must be validated. Keep context-specific components, and learn pooling weights on held-out lineages rather than assigning identities to A/B/C from control correlations.
3. **Diagnose the response head.** Use nested public-source folds to test ridge penalties spanning 0.01 to 10 and a small response-amplitude grid. Fit source-noise shrinkage inside each training fold. Compare known-target transfer and unseen-target prediction separately. Select using corrected expression error, fidelity, overlap, reach and discrimination together, not average MSE alone. Record response spread and predicted DE counts to detect collapse toward zero.
4. **Isolate count emission.** Hold response estimates fixed. Compare V4 emission with one-cell mean/pooled matching and a small dispersion grid, preserving integer counts, detection and realistic library distributions. Begin with neutral control-to-control tests, then known source perturbations. Measure the intended versus realized response, false DE calls, cell variance and between-target spread. Apply the actual official metric implementation/version when possible; keep proxy scores explicitly labeled.
5. **Evaluate a transition neural model only after the expanded baseline.** Compare small perturbation-conditioned neural transition heads on the expanded genetic-perturbation sources; a full STATE training/adaptation remains a separate, larger experiment. Compare it with transfer at matched data coverage and the same decoder. Add a learned residual only if it consistently improves entire-lineage holdouts. Larger observational encoders alone are insufficient evidence for adoption.

## Validation and submission decision

Use entire-lineage holdouts with unperturbed query controls permitted, plus donor/batch separation and target-family holdouts. Tune only inside training folds. H1's nine well-populated evaluation targets are useful but too narrow to be the sole selection set. Previously inspected H1 data remain development data; a new split does not make them a fresh blind benchmark.

Use the same template cells and random seeds for paired candidate comparisons. Evaluate a coverage-only change before combining new pooling and decoding. Report target-level uncertainty and disagreements among held-out contexts. Preserve V4 until a candidate shows consistent supported improvements; a local improvement cannot guarantee a hidden score improvement. Public descriptions cover only disclosed implementations, not all private high-ranking teams. Scores and historical ranks should not be compared as identical experiments unless panels and anchors match.

## V7 implementation update

The approved follow-up is implemented by root-level `prepare_v7_sources.py`, `audit_v7_source_offsets.py`, `fit_v7_response.py`, `train_v7_transition.py`, `evaluate_v7_response.py`, `fit_v7_submission.py`, `build_v7_submission.py` and `validate_v7_submission.py`. V7 uses a target/batch-stratified sample across the complete pinned public release, not every raw cell. The extraction found a gap in the converted expression table that invalidated later metadata offsets; it now audits fragments, excludes unavailable cells and corrects offsets, then verifies cell IDs and full libraries before projection. The expanded source extraction and all nine generated-cell comparisons are complete. The frozen calibration guards selected the V4 transfer fallback: none of the new estimators demonstrated an acceptable improvement. Upload readiness is recorded only after independent package validation passes. See [V7_MODEL_STRATEGY.md](V7_MODEL_STRATEGY.md) for the full decision and numerical results.

The four-penalty ridge study is implemented in `train_v7_ridge.py`, with nested lineage tuning and a conservatively masked residual comparison. A faithful Arc STATE reproduction is not implemented in this iteration; the neural comparison uses three small response architectures with inner-fold epoch selection. No claim about their hidden performance is made.

The V7 follow-up is now complete: official packaging and the independent full-payload audit passed for all 360,000 cells and 900 groups. [V7_SUBMISSION.md](V7_SUBMISSION.md) records the upload-ready file and comparisons. The selected effect arrays exactly equal V4; the package is a new sampling realization, not a demonstrated improvement. It is published as `k0pkYuJRcw3C0U8K0NgL`, overall **-0.051386751**, slightly below V4. See [V7_LEADERBOARD_RESULTS.md](V7_LEADERBOARD_RESULTS.md).

## Post-publication research and metric audit

The fresh review in [VCC2026_RESEARCH_AND_GAPS.md](VCC2026_RESEARCH_AND_GAPS.md) supersedes the next-step ordering above. A hash-verified retrospective audit of all 94 V7 profile caches demonstrated artificial discrimination from target-specific control origins: constant predictions score 0.854545-1.000000 historically versus 0.5 using a common control. Templates were also matched to actual treated batch proportions, unavailable for challenge predictions. Rebuild context/batch-compatible validation with target-independent templates before selecting another candidate. The expanded models also failed expression guards; this finding alone does not justify adopting them. Historical scores and packages are preserved. Reproduction and limitations: [TRANSFER_VALIDATION_AUDIT.md](TRANSFER_VALIDATION_AUDIT.md).


## V8/V9 corrected-validation implementation

The historical common-reference audit is now integrated into separate root pipelines. V8 compares V4 response amplitudes, pooled multi-source estimates and adaptive shrinkage with three decoders. It uses control-only technical strata, common reference origins, globally disjoint calibration/confirmation targets and target-independent templates. Whole-gene expression jackknifes and direction-yield budgets follow the pinned official source review in `analysis/v8_official_metric_review.json`. The previous local PDS values must not drive new choices.

V8's reduced-variance decoder was rejected by disjoint challenge-control tests: 968-2,403 false DE calls without perturbation. Full matching passed the neutral gate but lost supported-target fidelity on confirmation. V9 therefore isolates count initialization: start from V4 zero-aware simulated perturbations and calibrate only a selected fraction of population mean/pooled expression, preserving activated zeros and simulated libraries. It reserves further count-validation targets, with >=150 actual cells for HCT/HEK and >=60 for H1. Selection is frozen before confirmation and all outcomes are recorded without reselection.

[The V8 report](V8_SUBMISSION.md) and [V9 report](V9_SUBMISSION.md) preserve candidate-level metric tradeoffs, neutral checks, counts, hashes, reproduction and upload steps. An upload-ready package establishes format/integrity; public proxies cannot establish a higher hidden score. Broader matched interventional training, source platform/time/dose calibration and a successfully held-lineage validated response model remain research gaps. The control-expression reconstruction networks do not resolve causal transfer. V4 remains the best published candidate until a new leaderboard result demonstrates otherwise.

<!-- V9 published decision -->

## V9 publication update

V9 published at -0.070447281, below V4 (-0.050673910). Five scaled metrics improved, but fidelity lost -0.237327864 and outweighed their gains. The public mean-fidelity drop was already visible; the reference-supported weighted guard was insufficient. Retain the corrected scorer, preserve V4 as best, and require both overall and supported-target fidelity protection plus explicit metric availability before promoting another decoder. See [V9_LEADERBOARD_RESULTS.md](V9_LEADERBOARD_RESULTS.md).


## V10 population transport experiment

V10 makes a new response estimate from all five measured source families, with knockdown-dose correction and empirical Bayes shrinkage after pooling. It isolates real-template versus compact synthetic emission and a response-selective variant. The October 7 public Sidechain result shows a similar fidelity collapse when real control templates replaced synthetic cells, despite local gains; this motivates a new dispersion test rather than assuming more realistic variation improves the competition metric.

The original frozen V10 calibration snapshot (before the validation-depth adapter correction described below) selects `dose_compact`, with full target-mean fidelity 0.689284 versus V4/legacy 0.470546 and paired-finite fidelity 0.687225 versus 0.470546. Its expression ratio is 2.503636 versus 4.929231 on public source panels, and balanced discrimination is unchanged at 0.520833. DE overlap and reach decline. Confirmation must be assessed without reselection. These raw public results are not hidden leaderboard predictions.

This experiment explicitly changes the prior neutral-preservation policy: compact populations produce 2,039-3,265 neutral DE calls and fail the earlier safety threshold. That bias is measured and documented, not declared solved. The candidate is an experimental scoring-oriented population predictor, not a claim of faithful biological single-cell variance. Source-type and assay transfer still limit performance. See [V10_STRATEGY.md](V10_STRATEGY.md) and the generated [V10_SUBMISSION.md](V10_SUBMISSION.md).


### Final V10 depth correction and evidence

The validation adapter supplied full-population depths for a 400-cell template, incorrectly using the first 400 population depths. Original snapshots are retained as historical evidence; production already uses correctly indexed full-control depths. `audit_v10_validation_depths.py` repaired only the validation adapter and reassessed the frozen candidate without retuning or reselection. These are revisited public development cases, not fresh blind confirmation.

Corrected full independent reference comparisons (768-2,144 controls) give balanced panel fidelity 0.525383 -> 0.574049, discrimination 0.701389 -> 0.750000 and expression-error ratio 2.913312 -> 1.833134. Both full-target and paired-finite fidelity guards pass; paired six-metric proxy utility improves by 0.122538. Reach declines (0.611467 -> 0.385929), and fold-change error increases slightly (0.869425 -> 0.881051). The gain in fidelity is much smaller with full references than with 400 references. Neutral-population bias remains unresolved. This supports trying the candidate, not claiming that it beats V4 on the hidden leaderboard.

The complete experiment, packaging memory workaround, corrected diagnostics and upload instructions are documented in [V10_STRATEGY.md](V10_STRATEGY.md) and [V10_SUBMISSION.md](V10_SUBMISSION.md).

Full-reference fidelity also declines for HCT116 (0.664027 -> 0.587194), while H1 and HEK293T improve. Aggregate guard passage is not evidence that every unseen lineage will improve.

<!-- V10 published decision -->

## V10 publication update

V10 is published at **0.116754121**, beating V4 by **+0.167428032**, and is now the best project submission. Discrimination and fidelity account for 97.42% of the net gain; expression and overlap also improve, while reach and fold-change accuracy worsen. Preserve V10 and use controlled ablations to separate response pooling from emission. Earlier V4-best statements describe historical decisions. Neutral-control bias remains unresolved. See [V10_LEADERBOARD_RESULTS.md](V10_LEADERBOARD_RESULTS.md).
