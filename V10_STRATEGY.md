# V10: measured intervention responses and population transport

## Why later versions did not beat V4

V4 is still the best published result (-0.050673910). V9 scored -0.070447281 despite improving five scaled metrics: raw expression error fell from 1.247733 to 1.011257, but raw direction fidelity fell from 0.309235 to 0.238046. Fidelity contributed -0.039555 to the overall change, outweighing +0.019781 from the other five metrics. These are published comparisons, not estimates of V10.

V7/V8/V9 production responses ultimately retained V4's biological effects. Testing new neural architectures did not produce an accepted new response estimator. Observational prediction of other genes from the 300 inputs does not establish the causal effect of silencing an input. Our original discrimination validation also had target-specific control origins; the repaired tests use one shared control and return exactly 0.5 for constant predictions. V9's later guard protected only a supported subset of fidelity and missed its decrease in the full target mean.

## What public competitors disclose

[Sidechain's current model notes](https://saberhq.com/sidechain/models/) report a best score of 0.1626 with pooled source responses, shrinkage after pooling, knockdown-strength correction, and synthetic population emission. A new October 7 real-control-template variant reportedly fell to -0.0039 despite gains on local tests, with fidelity responsible for the loss. This is participant-reported evidence, not an independently reproduced result or an explanation of every leading team's methods.

[AtlasShift](https://github.com/kaipengm2/Virtual-Cell-Challenge-2026) reports 0.1545618019 and historical rank 82 on September 10. Its code separates average CPM from pooled profiles, combines public sources, and rounds integer counts subject to both moments. We reuse the already reviewed, hash-verified integer fitting routine, not the author's claimed performance.

[Arc's October 1 data article](https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge-2026) describes consistent Flex profiling and curated perturbations with at least 80% median knockdown. Public source screens differ in platform, context, depth and knockdown strength. The resulting transfer gap is a plausible limitation, not a measured explanation of our hidden-panel errors.

[VCLMU](https://arxiv.org/abs/2610.04475) represents interventions as transitions in latent mechanisms. It does not provide evidence that our small observational networks solve this competition. V10 uses actual measured genetic interventions and controlled emission; it is not a reproduction of that foundation model.

## New experiment

`v10_methods.py` compares six isolated strategies. V4/legacy remains the paired reference. V4/compact isolates emission. Pooled/compact changes the biological response estimator. Dose/compact also corrects measured source knockdown to 90%, with gains bounded to 0.8-1.8. Invariant/compact subtracts half of each source's common response before pooling to test sensitivity to shared stress/platform shifts. Dose/selective compresses variability only in proportion to response strength.

The new pooled response uses five source families (K562, CD4, H1, HCT116, HEK293T), existing measured-gene masks, weights [2, 0.25, 2, 1, 1], bounded confidence, and target-specific mixture-of-normal empirical Bayes shrinkage **after** pooling. Sources missing a gene abstain. Existing source estimates already include source-specific preprocessing; this is not a full raw-data reproduction of Sidechain or AtlasShift. The fit uses no challenge treated cells.

Compact emission starts from a new population profile with half-Poisson noise amplitude (one-quarter Poisson variance before fitting). It retains sampled observed cell depths and uses the reviewed AtlasShift integer transport to fit both average per-cell composition and pooled composition. Synthetic cells do not copy individual control expression vectors. Selective emission retains individual variation in genes with weak predicted response. Own-target reduction is handled separately and excluded from validation metrics.

## Validation and limits

`prepare_v10_validation.py` freezes control-defined technical strata from V8, new independent control draws, 400 predictions and 400 references per stratum/seed, and distinct actual source cells. Every response benchmark excludes the entire evaluated source family. Calibration and confirmation target sets are globally disjoint within this experiment. HCT116/HEK293T count targets exclude previously tested V7/V8/V9 targets. H1's small 25-target set is exhausted, so reused H1 targets are explicitly development evidence. Public source summaries were previously available: this is not pristine blind validation.

`evaluate_v10_submission.py` computes actual asymptotic Wilcoxon/BH calls, corrected expression profiles, shared-control discrimination, direction fidelity/reach, significance overlap and DE fold-change error. Own targets are excluded; discrimination excludes all 300 panel targets. Undefined metrics and denominators are retained and recorded. `select_v10_submission.py` compares a common set of complete six-metric panels and guards BOTH full target-mean and paired-finite fidelity within 0.01 of V4, supported-budget fidelity within 0.02, PDS within 0.03 and expression error within 10%. Equal family/context weighting reduces domination by larger sources. The calibration choice is frozen before confirmation; failed confirmation is reported without reselection.

**Neutral population bias is a material limitation.** Compact emission produced 2,039-3,265 significant DE calls against independent challenge controls even with zero response, with about 17% of original CPM variance. It fails the former control-preservation gate. Legacy and selective emission produced zero such calls in these six controls-only comparisons. V10 explicitly measures this bias and treats the compact model as an experimental competition population predictor; it must not be interpreted as a biologically faithful control reconstruction. V10 does not silently declare that earlier safety criterion passed.

Count panels contain only 2-4 targets and source-treated groups often have 150-200 cells; their discrimination and DE power differ from the 300-target, 400-cell hidden panel. Raw public utility is not the official anchored leaderboard score. The new artifact's integrity can be verified locally; improved leaderboard performance requires submission.

## Reproduction

Run from Windows Command Prompt in the project root:

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe run_v10_pipeline.py --workers 4
```

To preserve completed frozen stages and continue:

```bat
.\.competition_env\Scripts\python.exe run_v10_pipeline.py --resume --workers 4
```

The pipeline runs scientific invariants, prepares the frozen plan, evaluates calibration, records neutral bias, selects, evaluates confirmation once, fits all public sources, independently recomputes the model and checks splits, performs a 24-group preflight, generates all 900 groups, calls official `vcc prep`, compares every packaged count/index/pointer and label, and exactly replays 12 groups. Source preparation prerequisites and original input hashes remain recorded in V7/V8 and source provenance. Numerical caches and large submission artifacts are ignored by Git; code and small reproducibility records are tracked.

See `V10_SUBMISSION.md` after the full pipeline completes for the actual selected model, numerical comparisons, verified package hash and upload command.


A Git clone contains the frozen small records but excludes large numeric caches. Restore the competition inputs, hash-verified V4 source components and V7 source cache using the existing source preparation instructions first, then use `--resume`. `restore_v10_model.py` reconstructs a missing V10 model cache without changing its frozen manifest and requires its exact SHA-256 to match. It refuses to overwrite an existing cache. NumPy/runtime changes that alter compression or numeric output are reported, not silently accepted. Fresh large prediction files are regenerated and audited from this restored model; upload-package timestamps may change, while integer prediction replay remains the criterion for reproducibility.


## Validation depth-adapter correction

During the final audit we found that compact source-validation calls passed full-population depths together with a 400-cell template. `emit` indexed the population's first 400 depths, instead of the template's sampled depths. This affects the original calibration/confirmation evidence; those records and their frozen code are preserved as historical snapshots. Production calls use the full provided control matrix and its correctly indexed depths, so the generated submission is unaffected.

Root `audit_v10_validation_depths.py` supplies the intended template depths for compact emission, preserves legacy's population mean for its zero-aware rule, verifies exact generated depth equality on every corrected case, and reassesses the frozen candidate under both seeds with 400 and full independent reference populations. It does not retune or reselect. `--production-only` checks all 360,000 production depths against their recorded sampled controls after independent payload validation. The final report prioritizes these corrected results and explicitly labels original snapshots. The original `check_v10_reference_size.py` diagnostic also uses the historical adapter; its results are superseded by this corrected audit.


## Official packaging memory adaptation

The complete matrix has 2,459,676,692 nonzero counts, crossing SciPy's 32-bit-index threshold. The initial official prep could not allocate a further 18.3 GiB int64 index array on this Windows host. `package_v10_memory_safe.py` streams a local staging H5AD with float32 **unchanged integer-valued counts** and int64 CSR indices/pointers. Every staged value and coordinate is compared to the original before the unchanged official CLI runs with all normal checks. This avoids sparse index-widening and datatype-copy transients. The original count matrix and frozen generation code remain unchanged. The final independent audit compares every packaged value/index/pointer to that original matrix.

`generate_v10_matrix.py` calls the frozen generator as a function, validates the full matrix, and lets the pipeline use this safe packer directly on reproduction. The earlier packaging failure is retained in logs and official stderr; no shape, density, integer-count, target or cell-count check is disabled.


Corrected full-reference fidelity varies by held-out family: H1 0.447434 -> 0.599574, HCT116 0.664027 -> 0.587194, HEK293T 0.472062 -> 0.535381. The aggregate improvement therefore does not mean uniform transfer improvement. Larger references sharply reduce the apparent gain relative to sampled 400 controls. The generated report includes both reference sizes and the family breakdown.


`audit_v10_provenance.py` links the corrected score files, production-depth provenance, lossless staging record, independent payload audit and frozen code hashes. `document_v10_submission.py` requires those checks before writing the upload report or marking V10 ready in the README. No script uploads the package automatically.

The final audit records the exact Python, numerical-library and official CLI versions in `experiments/v10_population_transport/analysis/runtime_versions.json`. Use those versions when restoring frozen numerical caches.

<!-- V10 published result -->

## V10 publication update

V10 is published at **0.116754121**, beating V4 by **+0.167428032**, and is now the best project submission. Discrimination and fidelity account for 97.42% of the net gain; expression and overlap also improve, while reach and fold-change accuracy worsen. Preserve V10 and use controlled ablations to separate response pooling from emission. Earlier V4-best statements describe historical decisions. Neutral-control bias remains unresolved. See [V10_LEADERBOARD_RESULTS.md](V10_LEADERBOARD_RESULTS.md).
