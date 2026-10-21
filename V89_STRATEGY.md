# V89: release test of the earlier unsubmitted V86 projection

V89 is a release and repeatability test of V86's original `interval25` candidate. It does not claim a new learning method. V88's converged mixture fit improved source likelihood but failed prospective donor discrimination/overlap limits, so V88 is not packaged. Its complete two-seed negative results remain in V88_RESULTS.md and the original analysis directory.

The previously audited V86 candidate had improved five of six donor endpoints, with unchanged discrimination. V89 selected that known candidate before examining its new second-generation-seed result. The original V86 public predictions, model audits and seed0 count records are reused byte-for-byte and explicitly labelled reused. The donor and public studies have already informed project decisions; this is not fresh biological validation.

The fixed candidate retains V86's 100-update covariance-mixture fit, its normal-moment interval projection, quarter-strength bounded correction, and original V66 identity guard. These approximations are retained as the empirically more conservative configuration; fitting is not claimed converged, and posterior intervals are not claimed calibrated under biological domain shift. No numerical-optimization superiority is inferred from local biological accuracy.

The V89 protocol was frozen before the candidate's new second-seed test. Both seeds must satisfy the same bounded experimental-release limits used for V88: expression error within0.5%, PDS loss<=0.005, fidelity loss<=0.002, Jaccard loss<=0.0005, NMAE within1%, and at least two improved endpoints. Broad profile/PDS limits also apply. Original stronger criteria remain separately reported, including failed shuffled controls. Controls and other strengths cannot be promoted after failure. The exact thresholds, prior evidence and source hashes are in `experiments/v89_v86_projection/analysis/plan.json`.

Production, if permitted, uses all five measured source families and original V66 A/B/C effects, template identities, cell libraries and seeds. The V86 implementation is loaded in a separate Python module namespace with a V89 output folder; original source files and records are never modified. Float64 Cholesky arithmetic is substituted only after actual-source parity checks against generic Gaussian solves. Since the immutable inference computes controls alongside the primary, production real and shuffled priors are fitted, but only `interval25` is deployed. The unchanged V66 residual10 emitter produces the cells.

## Reproduction

Use the documented scientific Python3.12 runtime for scientific stages, and `.competition_env\Scripts\python.exe` for official packaging. Run from the root:

```text
v89_projection_release.py freeze
v89_projection_release.py adopt
run_v89_counts.py counts --offset 620000
mean_v89_release.py
run_v89_counts.py assess --offset 0
run_v89_counts.py assess --offset 620000
audit_v89_counts.py --offset 620000
decide_v89_release.py
```

Count jobs can use disjoint `--chunk K --chunks N` partitions. Seed0's original V86 audit is reused and labelled with zero current emissions. `finish_v89_release.py` waits for the four disjoint count workers used in this run, applies the fixed decision, and executes production only if allowed. For a fresh reproduction use a separate output checkout/folder; existing files deliberately cannot be overwritten.

The production plan hashes its adapters before fitting. `fit_v89_submission.py freeze` then the all-source real/shuffled finite-iteration prior, `audit_v89_production.py`, and `fit_v89_submission.py fit` prepare responses. The driver then runs preflight, all900 groups, official prep and independent validation. `build_v89_submission.py` streams360,000cells across18,533genes with exact original libraries; `package_v89_submission.py` stages integer-valued counts losslessly and invokes official prep; `validate_v89_submission.py` compares the complete decompressed payload, every original template/library, and12 fully regenerated400-cell groups.

The output, if all integrity checks pass, is `experiments/v89_v86_projection/prediction.vcc`. A passing file format is not proof of biological accuracy. V89 remains explicitly experimental unless the separately reported stronger criteria pass. No script uploads or pushes to GitHub.
