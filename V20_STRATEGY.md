# V20: paired-control Stack response ensemble

## Purpose and fixed design

The pretrained Stack population model is used as a bounded response expert alongside the V16 predictor. It is not asked to replace all genes or to emit the submission directly. The development experiment compares 196 complete integer-count populations: seven arms, 14 context/target pairs and two fixed seeds. The independent algorithm reserve contains eight other context/target pairs, with two seeds each. Public biology has partly been used in earlier work; the reserve is an algorithm check, not a new blind biological dataset.

Author code is pinned to ArcInstitute/stack commit `cacc2e4b09435c3e536d46237d10b50f222dd144`. Stack-Large-Aligned weights are pinned to model revision `b09f085dac03d170b078a5c72f550ae93686e544`, checkpoint SHA-256 `f93cf6f42f36c8a85dc570d92e801c1fc3e1f45d55741e6bb250d178a6b6ad36`. The actual configuration has approximately 217 million parameters and 15,012 output genes. Of these, 13,078 occur in the competition axis. Unsupported coordinates retain V16 response multipliers. Follow the archived code/model license and model policies; these are not unrestricted commercial weights.

The whole query family is excluded from response fitting and source prompts. Its controls may be supplied as the query. Stack pretraining overlap with public studies has NOT been independently excluded; prompt exclusion does not establish unseen-pretraining-family generalization. No observed query treatment is an inference input.

## Inputs and model execution

For each target, select a source family with at least 40 treated cells, using direct pairwise-valid Spearman correlation of source and query control mean CPM. Exclude all 300 panel genes from this comparison; do not add a pseudocount or rank log-transformed values. Pair each selected source treatment cell with a source control from the same batch. Repeated conditioning controls are disclosed, rather than counted as independent treatments.

The shuffled-target negative control uses the same source family and source-size bin. Real and shuffled prompts use the same prefix size, `min(real_count, shuffled_count, 400)`. Query input consists of 400 fixed genuine control cells. Missing vocabulary coordinates are zero-padded with explicit measured masks; these zeros do not imply biological absence.

Run strict full checkpoint loading, FP32 and batch size four. The one-step mean adapter exactly matches the author's `get_prediction` conditional mean and dispersion arrays, bitwise, on all 400 by 15,012 smoke-test coordinates; repeat inference also matches bitwise. It removes negative-binomial sampling and multi-step MDM iterations. This is a bounded one-step expert, not a reproduction of the author's five-step demonstration. Native sampling smoke checks were also performed. See `mean_adapter_parity.json` and `smoke*.json`.

## Response correction

Generate both source-control-conditioned and source-treatment-conditioned query means. Normalize each cell within their common measured output vocabulary. Let `L` be log2(mean treatment / mean generated control), with epsilon 1e-9. Variance combines conditional-mean cell variation with the analytic NB variance `mu + mu^2/theta`, using an independent-arm approximation. This is a shrinkage diagnostic, not an FDR estimate.

```
w = L^2 / (L^2 + variance + .0025)
    * .04 / (.04 + variance)
    * n / (n + 100)
correction = .25 * w * clip(L - V16_effect, -.5, .5)
```

Apply the existing query-gene and own-target expression gates to the correction only; do not gate the baseline twice. Preserve the own-target effect and the V16 fallback. The added correction is at most .125 log2 before gating. A global composition normalization means unchanged fallback multipliers do not guarantee unchanged emitted gene counts.

The selected mean arm uses the unchanged V16 residual10 integer emitter: .65 Poisson, .25 mean and .10 empirical-control components, followed by dual-moment fitting. Every emitted cell has exactly its original control library size. Mean and bulk probability targets are fitted approximately and their residual L1 errors are recorded; the frozen plan's phrase 'exact mean/bulk count constraints' is imprecise and must not be read as an exact achieved marginal equality.

## Prospective assessment

Guards are frozen in `analysis/plan.json`, SHA-256 `afdb437428290453b714b613ad97ff965a16ddd4cf14d59d37302ee2ef6d911e`. Assess full emitted counts against independent controls. PDS excludes all 300 panel genes; DE excludes the own target. Require the available overall and worst-family guards, at least three improved metrics including DE, and improvement over amplitude and shuffled controls on at least two of fidelity/reach/Jaccard. NMAE is missing where fewer than ten true DE genes are available; it is never replaced with zero. Missing metrics prevent a complete six-metric confirmation claim.

Only a development-selected candidate may enter the once-only reserve. A reserve failure must not be followed by parameter adjustment on that reserve. These raw public metrics do not reproduce the hidden competition anchor-scaled score or guarantee ranking.

## Reproduction

Use `.competition_env` for the unchanged count-scoring numerical stack. Use the isolated Python 3.12 Stack runtime recorded in `analysis/inference_runtime_freeze.txt` for checkpoint inference. Run `setup_v20_stack_runtime.py --install` to create a standalone local Python 3.12 runtime with the exact recorded package pins; `--verify-runtime <python.exe>` checks an existing runtime. Runtime and asset setup are separate from model selection. Downloaded weights, raw data and numeric caches are ignored by Git; their hashes are recorded.

```powershell
.\.competition_env\Scripts\python.exe collect_v20_stack_assets.py
.\.competition_env\Scripts\python.exe prepare_v20_inputs.py development
# Run the strict checkpoint/mean-adapter parity smoke test in the Stack runtime.
.\.competition_env\Scripts\python.exe prepare_v20_experiment.py
# Stack-runtime command, once per family H1, HCT116, HEK293T:
python infer_v20_stack.py development --family H1
.\.competition_env\Scripts\python.exe run_v20_experiment.py development --workers 6
.\.competition_env\Scripts\python.exe audit_v20_inputs.py
.\.competition_env\Scripts\python.exe audit_v20_experiment.py development
.\.competition_env\Scripts\python.exe assess_v20_experiment.py development
# Only after a passing selection:
.\.competition_env\Scripts\python.exe prepare_v20_confirmation_inputs.py
python infer_v20_stack.py confirmation --family H1
# Repeat inference for HCT116 and HEK293T before scoring.
.\.competition_env\Scripts\python.exe run_v20_experiment.py confirmation --workers 6
.\.competition_env\Scripts\python.exe audit_v20_experiment.py confirmation
.\.competition_env\Scripts\python.exe assess_v20_experiment.py confirmation
```

Existing frozen plans/assessments deliberately refuse overwrite. Preserve them and their hashes. See V20_RESULTS.md for outcomes and the release decision.

The original repeated input preparer compares Python tuples with JSON lists and stops on an equality assertion. The separate reserved-input helper reuses the frozen choices, exactly reconstructs all same-batch source rows and count values, and records its own manifest. It changes no numerical model code or frozen plan.

## Production after reserved approval

`v20_production.py` requires the passing mean-only development/confirmation and independent source/case audits. It checks its source-choice implementation against all 900 frozen public target/family choices before preparing A/B/C. All three public source families may be used for these new contexts, using the same rank comparison, minimum cell count, source-size bins and matched-prefix size. The shuffle rule still fixes the validated prompt size, but shuffled predictions are not generated in production. A 400-cell guide/depth-stratified query control sample is fixed with seed sequence `[20261020, context_index]`. The global mean/bulk inputs still use all 18,400 original controls per context.

Production uses only the validated Stack mean fields and unchanged V16 residual10 emitter; V21 active variation is excluded. Source rows and same-batch controls are recorded. A fast exact sparse vocabulary aligner is checked against the frozen original on each private query. `production_plan.json` freezes code, input hashes and the seed rule before inference. `infer_v20_submission.py` strictly loads the checkpoint once, infers 900 target/context experts and checkpoints resumable numeric arrays every 20 targets. `fit_v20_submission.py` applies the exact selected correction to frozen V16 arrays, bounds changes by .125 log2 and records all changes/fallbacks.

```powershell
.\.competition_env\Scripts\python.exe v20_production.py
# Use the isolated Stack runtime for this command:
python run_v20_submission_inference.py
.\.competition_env\Scripts\python.exe fit_v20_submission.py
.\.competition_env\Scripts\python.exe build_v20_submission.py --preflight --workers 4
.\.competition_env\Scripts\python.exe build_v20_submission.py --workers 4
.\.competition_env\Scripts\python.exe validate_v20_submission.py
```

The streaming builder retains original V16 control sampling/count seeds and may reuse only independently verified V16 groups whose effects are bitwise unchanged. Packaging losslessly widens sparse indices and stores exact integer-valued float32 counts, then runs the unchanged official VCC prep. This avoids large CSR conversion copies above 2^31 nonzeros. The validator independently reconstructs all 900 correction vectors, compares every packaged count/index/pointer/axis, verifies all original control identities/guides/library sizes and exactly replays 12 groups, including the largest response changes. A package is declared ready only after all these checks finish. See V20_SUBMISSION.md for the final artifact receipt and upload command.

The full per-source production prompt matrices and repeated per-prompt row manifests are regenerable caches. They remain local and are SHA-256 referenced by `production_plan.json`; the compact private source-choice records, model metadata, audit receipts, group summaries and code are versioned. The 360,000-row prediction control provenance is also local, with its hash and full audit retained in Git. Original raw files remain required to reproduce the model.

A production metadata check caught an object-dtype export of gene/target axes before any expert output. The failed production manifest and original exporter files are preserved under `production_plan_before_axis_dtype_fix.json` and `production_axis_export_snapshot/`. `repair_v20_axis_export.py` reads official CSV string axes, never unpickles an array, and verifies every other numeric array and all source files remain unchanged. `production_axis_export_repair.json` records the repair; the current production manifest includes its hash. This is a safe export correction, not a model/parameter adjustment after reserved scoring. Fresh production code writes plain-string axes directly and does not require the one-time repair.

Run the additional integrity gates before declaring an artifact ready:

```powershell
.\.competition_env\Scripts\python.exe audit_v20_production_inputs.py
# After build_v20_submission.py and validate_v20_submission.py:
.\.competition_env\Scripts\python.exe release_v20_submission.py
```

The Windows Python 3.12 inference runtime hit a native pandas/Arrow string-construction access violation while reading the specification CSVs after loading Torch. `run_v20_submission_inference.py` applies `pd.options.future.infer_string=False`, the same metadata flag used in the author-output parity smoke test, before calling the frozen inference function. Axes are explicitly converted to Unicode; no scientific code, weights, numeric inputs or sampling parameters change. The diagnosed stack trace is retained locally, and the wrapper hash/flag/reason are recorded in `production_runtime_compatibility.json`. Use this wrapper for production inference.

For the remaining build stages, `.\.competition_env\Scripts\python.exe finish_v20_submission.py --workers 4` checks completed experts and the raw-input audit, then runs fitting, preflight, generation, packaging, full validation, release and documentation in order. `--wait-inference` permits preparing this runner while inference finishes; it starts no fit until all three 300-target checkpoints are complete. Verified completed stages may be reused when resuming after a failure; inspect the error before restarting.
