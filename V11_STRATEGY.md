# V11: consensus-supported coverage and controlled dispersion

V10 is the published baseline: overall 0.116754121 on `vcc2026-val-1`, entry `CYpAqJP7bUYqjA2ouD7U`. Its discrimination and fidelity gains beat V4, but raw fold-change error and reach worsened. V10's model, package, numerical code and validation snapshots are preserved.

## Research and design

The fresh [Sidechain model review](https://saberhq.com/sidechain/models/) discloses separate experiments for source shrinkage, synthetic-cell noise and mean versus pooled response strength. Its local gains sometimes failed on the leaderboard. No newly disclosed recipe beyond SER15 was found in this review. These participant reports motivate ablations, not a claim to reproduce their score. V11 reuses the pinned [Arc metric review](https://github.com/ArcInstitute/cell-eval2/tree/0ce25cf4c495e76e60dfca4ddc8f822507aa0e18) and hash-verified [AtlasShift integer fitter](https://github.com/kaipengm2/Virtual-Cell-Challenge-2026).

V11 freezes nine comparisons before evaluation:

| Candidate | Response | Noise amplitude | Pooled response multiplier |
| --- | --- | --- | --- |
| v10 | Exact V10 dose-corrected empirical Bayes response | 0.50 | 1.00 |
| v4_compact | Historical V4 transfer with compact emission, diagnostic only | 0.50 | 1.00 |
| dose125 | V10 response multiplied by 1.25 | 0.50 | 1.00 |
| consensus20 | Conservative consensus rescue | 0.50 | 1.00 |
| consensus40 | Stronger consensus rescue | 0.50 | 1.00 |
| noise35 | V10 response | 0.35 | 1.00 |
| noise65 | V10 response | 0.65 | 1.00 |
| consensus20_noise65 | Conservative rescue and higher noise | 0.65 | 1.00 |
| bulk085 | V10 per-cell response; weaker pooled response | 0.50 | 0.85 |

`v11_methods.py` starts from V10's dose-calibrated family pool and its empirical Bayes estimate. For each target/output gene, it calculates source-weighted sign agreement, requires at least two measured source families, and caps confidence using the pooled signal-to-noise ratio. It adds 20% or 40% of the gap between the pooled estimate and its shrunken value, multiplied by that support gate. Single-family, conflicting or uncertain responses are not rescued. Missing measurements abstain. Rescue preserves the pooled response sign; it does not use the held lineage's treated values.

The emitter retains synthetic populations rather than copying individual controls. Mean CPM and pooled composition are fitted separately with integer transport; every predicted library exactly matches a sampled control depth. Noise changes the mixture of expected and Poisson counts. Own-target reduction remains 90%, separate from downstream responses, and the own gene is excluded from metrics. A wrong full-population depth vector for a sampled validation template is rejected explicitly.

## Validation

`prepare_v11_validation.py` freezes control-only strata from the previous experiments. Original V10 calibration targets are explicitly reused development. V11 confirmation uses previously untested HCT116/HEK293T count targets, excluding all V7/V8/V9/V10 count-tested targets. H1 has only 25 targets and remains reused development biology. Public source summaries have been inspected, so this is not pristine blind validation. Development and confirmation target names are globally disjoint within V11.

Every source test excludes the entire evaluated family from pooling and shrinkage. Templates are target-independent, references are independent controls, and actual treated cells are distinct. There are 400 predictions, 60-200 observed H1 cells or 150-200 HCT/HEK cells, two paired seeds, and both 400 and full independent reference populations. The baseline V10 estimator and emitter are checked for exact numerical equivalence in focused tests; tests also cover source exclusion, support requirements, raw integer counts, deterministic replay and depth-axis correctness.

`evaluate_v11_submission.py` uses actual asymptotic Wilcoxon/BH calls, common-control discrimination, direction fidelity/reach, DE overlap and fold-change error, plus signed expression-error corrections with explicit missingness. Constant predictions must score chance discrimination. Full references receive 70% of the raw public utility and sampled 400 references 30%; families and technical contexts are equally weighted. This utility is not the official anchored score.

`select_v11_submission.py` protects V10 rather than V4. On both reference sizes, full and paired-finite fidelity may fall at most 0.01, supported fidelity at most 0.02, PDS at most 0.03, reach at most 0.03, and NMAE may increase at most 0.015; expression error may increase at most 10%. Full-reference fidelity for each source family may fall at most 0.03. Selection requires all three families and identical complete six-metric panels. The V4 compact ablation is not eligible. The winner is frozen before confirmation; a failed confirmation stops model fitting and is not retuned into a success.

The metadata audit found only two planned H1 groups with at least 400 distinct treated cells; HCT/HEK groups contain 150-247. Before selecting a model, `audit_v11_cell_power.py` therefore adds a separate 400-cell DE guard: development on ANKZF1 and frozen-candidate confirmation on SLIRP, each with both seeds and reference sizes. Fidelity may fall at most 0.02 and NMAE may increase at most 0.02 relative to V10. This amendment was based on metadata availability before inspecting the 400-cell results. It cannot resolve the power mismatch for the other lineages; no cells are duplicated to inflate significance.

## Distribution bias

`benchmark_v11_neutral.py` measures zero-response populations against disjoint supplied A/B/C controls, with no own-target knockdown, two seeds, and 400/2,000 references. Compact populations remain biased representations of biological variance. The noise-0.65 diagnostic reduces neutral false DE calls relative to V10; noise-0.35 increases them. Neutral results are reported rather than falsely claimed to pass the former control-preservation gate. A higher leaderboard score would not establish faithful biological single-cell variance.

## Reproduction and readiness

All new scripts are in the project root. Input/cache preparation remains documented in earlier V7/V10 reports. Large source caches and generated matrices are excluded from Git.

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe run_v11_pipeline.py --resume --workers 4
```

The pipeline checks scientific invariants, frozen development and confirmation decisions, all-source fitting, model provenance, preflight groups, complete generation, unchanged official VCC packaging and independent full-payload auditing. It never uploads. `restore_v11_model.py` recreates an ignored model cache only if it exactly matches its frozen SHA-256. `package_v11_memory_safe.py` losslessly stages integer-valued float32 counts and int64 indices to avoid conversion copies above 2^31 nonzeros; all counts and coordinates are compared before the unchanged official CLI runs.

`validate_v11_submission.py` independently recomputes production effects, checks source exclusions and split provenance, compares every packaged count/index/pointer/axis, verifies every cell's original control ID/guide/library and replays 12 complete target/context groups. `document_v11_submission.py` requires those checks before announcing readiness. The generated `V11_SUBMISSION.md` records the actual selection, comparisons, package hash, limitations and Command Prompt upload command. Improved hidden performance is unknown until submission.
