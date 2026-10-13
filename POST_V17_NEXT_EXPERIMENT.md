# Proposed next experiment after V17

Status: **research proposal, not an executed model comparison**. Prepared 2026-10-13. This does not modify V17's failed selection or any previous release policy. Freeze a new implementation, inputs, identities and thresholds before observing new model outcomes.

## First pilot: PIE as an independent response expert

**Question:** can a pretrained biological-knowledge model improve the identity and direction of DE genes beyond V16, rather than merely change overall response amplitude?

Use the pinned `PIE_replogle_wdataset` **K562-held-out** checkpoint (`055c7a2cabdaff5659121c64fa2b609c231fc5be`). It trained on HepG2/Jurkat/RPE1 and none of our 300 target identifiers. Do not use the cross-dataset checkpoint to claim H1/Orion held-family generalization: those sources are in its training set. Full provenance and code inspection are in [POST_V17_RESEARCH.md](POST_V17_RESEARCH.md).

### Phase 0: a bounded inference smoke test

1. Create a separate compatible Python 3.12 environment. Preserve `.competition_env`, current scorers and numerical dependencies.
2. Inventory required assets before downloading. The checkpoint is 949.6 MB, but a configured source array is 43.2 GB; establish a working source-mask or bounded-asset strategy first. Retain the checkpoint's original training split and evidence provenance. Official inference still needs its training-evidence assets.
3. Build one known example and one new controls-only context with no treatment labels. Check finite outputs, exact gene order, determinism, missing-feature masks and output units. `p_de` must never be substituted for a BH q-value.
4. Use a supported precision on the RTX 2070 SUPER (8 GiB; no native BF16). Start with CPU FP32 or a verified small GPU batch. Measure peak memory; do not assume the training configuration fits.
5. First restrict new effects to the 6,290 competition genes in the checkpoint's measured training vocabulary. Keep V16 for all remaining coordinates. The 18,262 available gene descriptions make a later output-vocabulary extension possible to investigate, not already validated.
6. Unknown context identity remains unknown. Supply real control expression and supported missing context-text inputs rather than guessed cell-line labels. Do not generate new text embeddings through a paid API for this pilot; published feature assets already exist.

Stop the pilot on incompatible axes/units, unavailable mandatory assets or unsupported source masking. Record the exact limitation; do not silently replace the model with an untested approximation.

### Phase 1: isolate the response expert

Use the existing verified H1/HCT116/HEK293T raw matrices and common-reference scorer. Keep V16 count generation, own-target handling, control sampling, count seeds and fallbacks fixed. A new public-trained expert is being evaluated here; no A/B/C treated truth exists.

Proposed fixed arms, to be frozen before scores:

| Arm | Purpose |
|---|---|
| V16 | Exact comparator |
| 75% V16 + 25% PIE LFC, on supported coordinates only | Test new effect directions/magnitudes |
| Same conservative blend with PIE DE-confidence weighting | Test whether its DE head adds information |
| V16 with target-wise L1 amplitude matched to the confidence-weighted blend | Separate learned gene selection from generic shrinkage |
| Target-shuffled PIE effects/confidence with the same support masks and blend | Negative control for target-specific information |

DE-confidence weighting initially uses a frozen bounded rule, not an interpretation of probabilities as calibrated FDR. Any learned probability calibration must use inner source-family folds only and match the cell-count/reference setting. Do not use the reserved outcomes to choose a threshold.

Measure full-axis raw PDS, corrected expression numerator/denominator, NMAE eligibility, fidelity, reach and Jaccard. Report DE call counts, sign errors, expert error correlation and genes on which PIE corrects or breaks V16. Evaluate complete emitted counts, not just the classifier on covered genes. Keep lower-is-better errors distinct from the official higher-is-better scaled scores.

Retain the existing available-metric/worst-family guards, require broad improvement including DE, and require the confidence arm to beat both shuffled and amplitude-matched controls. Missing NMAE remains missing. If development fails, stop rather than open reserved verification. The still-unopened V17 verification identities could be reused for a new once-only algorithm check after explicitly freezing their role; prior public biology is not blind independent biology.

The unchanged V16 emitter implies unchanged zero-effect count behavior. Do not claim that this first response-only pilot fixes the generator's null-DE or covariance problems.

## Second pilot: Stack with generated-control calibration

**Question:** does a pretrained conditional population model improve treated-state structure beyond V16 and beyond raw generative bias?

Use observed source treatment cells and matched source NTCs from a different evaluation family as prompts; query input contains only the held-family controls. The public drug tutorial motivates this setup, but does not establish CRISPRi success.

For each case, generate a treatment and a synthetic control with the same query cells and paired randomization. Compare their residual effect, then map that effect onto the real query-control baseline. Do not upload raw generated treatment counts on the assumption that generated-control cancellation also holds against real controls. First validate the residual-to-count transformation on disjoint NTCs and known public treatments.

The model covers 13,078 competition genes. Preserve a transparent fallback for the other 5,455 genes, account for differing library axes, and avoid interpreting zero-padded coordinates as biological zero predictions. Inspect public-pretraining overlap before labeling the result an unseen-family benchmark.

Proposed comparisons: V16; raw Stack as a bias diagnostic; Stack with generated-control calibration; matched shuffled source-treatment prompts; and a simple control-state-reweighting baseline. Measure response accuracy and joint distributions separately, with all six available count metrics, detection rates, gene covariance and frozen-control-space distances. The raw Stack arm is not automatically an eligible submission candidate.

Only after calibration and source-held-out gains survive a frozen verification should a separate V18-style artifact be built and fully audited. No model package is created by this proposal.

## CPU fallback: covariance-conditioned response forces

If checkpoint inference is impractical, test a small covariance-response transfer model using the existing raw caches:

- Fit one shared 8–32-dimensional program basis from training controls; preserve measured-gene masks and remove batch/depth nuisance effects.
- Estimate each target's force from source treatment mean shifts and shrunken source control covariances. Query response uses query control covariance, not query treatment data.
- Compare conserved transfer, real covariance, diagonal covariance and shuffled covariance/program labels. Keep V16 count generation initially.
- Test empirical state reweighting only after a control-support diagnostic. Occupancy changes are not identified absolute birth/death rates, and reweighting cannot create a new state outside control support.

CIPHER's treatment-informed fitted amplitude is an oracle diagnostic, never a deployable query input. This adaptation must earn its own held-family result.
