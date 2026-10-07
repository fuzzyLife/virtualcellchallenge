# Performance improvement plan

Research checked: 2026-10-07. This is a proposed experimental sequence, not a newly validated model or submission. Existing packages are unchanged.

## What our evidence says

V4 remains our best published candidate: overall -0.050674, versus V5 -0.113637 and V6 -0.056790 on the recorded validation panel. V4's scaled fidelity is -0.683855 and significant-set overlap is -0.029878. Its scaled expression score of zero is a clamped score, not zero error. These aggregate results cannot identify which context or target is responsible. Fidelity yield is not simply sign precision.

The control-expression neural networks learn observational relationships. A good prediction of other genes from the 300 measured genes does not establish their response to CRISPRi. The later embedding pilot does train on interventions: 889 source contrasts across 284 targets. However, it uses a single ridge penalty of 10, fixed response-label shrinkage and a fixed 25% blend. In its H1 public-source generated-cell test, Tahoe/rank reduces fidelity from 0.145431 to 0.057594, overlap from 0.014226 to 0.008027 and worsens the corrected expression ratio from 1.276690 to 1.283540. These are raw public development metrics, not leaderboard scores. Very weak predicted responses make underfitting a plausible explanation to test; this does not prove that lowering regularization will improve transfer.

Evidence and reproduction: [embedding pilot](EMBEDDING_RESPONSE_PILOT.md), [upstream review](UPSTREAM_REPOSITORY_REVIEW.md), [V6 results](V6_LEADERBOARD_RESULTS.md).

## Public approaches and material gaps

### Broader measured CRISPRi coverage

[AtlasShift](https://github.com/kaipengm2/Virtual-Cell-Challenge-2026) pools full prepared K562, HCT116, HEK293T, H1 and CD4 statistics, then reanchors responses on target-context controls. Its historical September 10 score is author-reported, not a current ranking or independently replicated result. Its preparation requires substantial downloads/scans.

Our K562 pseudobulk aligns only 7,681 of the 18,533 required output genes, despite covering 260 target genes. Our HCT116 and HEK293T data are small pilot batches, not comparable full intervention panels. Our updated H1 cache covers 25 targets. These are coverage gaps, not a reason to fill unavailable effects with zero or assume A/B/C identities from marker similarity.

### Response pooling and uncertainty calibration

[Sidechain's model log](https://saberhq.com/sidechain/models/) describes combining effects before adaptive shrinkage, source-confidence floors, measured knockdown calibration, and controlled count-emission experiments. It reports SER-14 at 0.1626 on October 3; its Tahoe-neighbor variant changed 0.1131 to 0.1132. These are participant-reported results and do not establish causality or guarantee gains here.

Our CD4 estimator already shrinks effects using standard errors, and K562 uses a noise threshold. The missing experiment is learning the shrinkage scale for a target's combined response and accounting for disagreement among sources, rather than applying only fixed source rules. Count-generation calibration is also not new to us: V5 tried separate mean-CPM and pooled-profile matching, alongside other changes, and scored worse. It therefore needs an isolated test, not automatic reuse.

### Perturbation-conditioned neural models

[Arc STATE](https://github.com/ArcInstitute/state) distinguishes a cell embedding model from a transition model trained to predict interventions. Its repository provides Replogle-Nadig genetic-perturbation training and whole-context zero-shot split examples. This is a relevant neural baseline, not proof of a superior current private leaderboard entry. Any comparison needs genetic-perturbation training, complete output-axis handling and conversion to the competition's raw count contract.

The [Stage 1 Lab](https://github.com/danielchen26/vcc2026-stage1-lab) explores rank-test detectability and count construction, but its README marks the biological response estimator unfinished. Treat its scorer claims as hypotheses to verify against all six official metrics, not as a demonstrated winning solution.

## Controlled sequence

1. **Expand source coverage first.** Extract all available panel-target and matched control cells from HCT116 and HEK293T, and replace the narrow K562 expression axis with full measured statistics. Record measured-gene masks, targets, guides, batches, donors, knockdown depths and uncertainty. Deduplicate controls across public splits. Freeze versions/hashes. Cache summaries so fitting does not require repeated full atlas reads.
2. **Create stronger transfer baselines.** Compare the original V4 estimator with full-source pooling, confidence floors and target-adaptive shrinkage. Test normalization by measured knockdown strength with conservative caps; this assumes transferable dose dependence and must be validated. Keep context-specific components, and learn pooling weights on held-out lineages rather than assigning identities to A/B/C from control correlations.
3. **Diagnose the response head.** Use nested public-source folds to test ridge penalties spanning 0.01 to 10 and a small response-amplitude grid. Fit source-noise shrinkage inside each training fold. Compare known-target transfer and unseen-target prediction separately. Select using corrected expression error, fidelity, overlap, reach and discrimination together, not average MSE alone. Record response spread and predicted DE counts to detect collapse toward zero.
4. **Isolate count emission.** Hold response estimates fixed. Compare V4 emission with one-cell mean/pooled matching and a small dispersion grid, preserving integer counts, detection and realistic library distributions. Begin with neutral control-to-control tests, then known source perturbations. Measure the intended versus realized response, false DE calls, cell variance and between-target spread. Apply the actual official metric implementation/version when possible; keep proxy scores explicitly labeled.
5. **Evaluate a transition neural model only after the expanded baseline.** Train or adapt STATE on the expanded genetic-perturbation sources. Compare it with transfer at matched data coverage and the same decoder. Add a learned residual only if it consistently improves entire-lineage holdouts. Larger observational encoders alone are insufficient evidence for adoption.

## Validation and submission decision

Use entire-lineage holdouts with unperturbed query controls permitted, plus donor/batch separation and target-family holdouts. Tune only inside training folds. H1's nine well-populated evaluation targets are useful but too narrow to be the sole selection set. Previously inspected H1 data remain development data; a new split does not make them a fresh blind benchmark.

Use the same template cells and random seeds for paired candidate comparisons. Evaluate a coverage-only change before combining new pooling and decoding. Report target-level uncertainty and disagreements among held-out contexts. Preserve V4 until a candidate shows consistent supported improvements; a local improvement cannot guarantee a hidden score improvement. Public descriptions cover only disclosed implementations, not all private high-ranking teams. Scores and historical ranks should not be compared as identical experiments unless panels and anchors match.
