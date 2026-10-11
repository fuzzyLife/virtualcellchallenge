# V14: conserved responses plus context interactions

This is a separately versioned research experiment. V13 remains the published reference. Numerical methods and target/control identities are frozen before outcomes in `experiments/v14_context_interactions/analysis/plan.json` and `execution_freeze.json`.

## Model

`v14_context_methods.py` retains the exact held-family V13 unknown-knockdown response as its conserved component. It fits two-dimensional context features from average-tied ranks of control expression on genes measured in every training family, excluding all 300 target-panel genes. Missing values are omitted from query ranks; insufficient coverage causes abstention. The feature basis, centers, scales and extrapolation radius use training families only. CD4 conditions are averaged within their family and never counted as independent contexts.

For each target and output gene, center dose-calibrated training-family responses and context coordinates using bounded measurement-confidence weights. Fit a weighted ridge slope, with at least three measured families. Compare penalties 2 and 8. The conditional correction is one-half of a clipped +/-0.5-log2 slope prediction, so its maximum contribution is +/-0.25. Add it to V13, then apply the existing query-expression gates and own-target exclusion. This can change response directions, unlike expression gates alone.

A fixed context-label permutation fits an equally dimensional negative control. A 1.1x V13 response comparator tests whether simple response amplification explains a gain. Zero, generic-response and shuffled-target summary controls distinguish target-specific signal. The shuffled-context control uses penalty 2; if penalty 8 advances, it matches dimensionality and the correction cap but not the penalty. It is explicitly not an exact penalty-matched ablation in that case.

This small model is not a reproduction of STATE, CellOT or a foundation model. The five public source families limit how well context interactions can be identified, and assay/context differences remain confounded. It tests known perturbations transferred to an unseen family, not unseen-target extrapolation. Training-family measurements of reserved targets are legitimate inputs for that task; reserved targets are excluded from summary evaluation and count-based model selection.

## Response-versus-population diagnostic

Select three targets per family using count availability and target-name hashes, before any V14 outcomes. Each has at least 400 genuine selected treatment cells. For each, randomly select 400 distinct cells and split into 200 treatment-calibration and 200 disjoint treatment-truth cells. Compare:

1. The held-family V13 response with the existing compact emitter.
2. Compact emission constrained to the measured mean and pooled proportions of the 200 calibration treatments.
3. Those actual 200 calibration treatments as an independent real replicate.

All three diagnostic arms use the same known treatment-calibration library depths. This deliberately holds depth fixed to study response/distribution limitations; even the first arm's depth input is diagnostic and differs from deployment. Neither treatment-derived oracle is eligible for selection or submission. No treatment cell is duplicated to claim extra independent observations. The real replicate is noisy, not a perfect upper bound.

Primary scoring uses full independent controls, the pinned CPU Wilcoxon/BH implementation, own-target exclusion, common-origin discrimination and jackknife-corrected expression. Oracle profile panels contain three targets; no private anchor scaling is available. A 200-cell diagnostic is not a 400-cell competition replica. Comparison with actual replicates helps expose the residual sampling/technical floor.

## Actual-count validation

Development contains three targets in H1, HCT116 and HEK293T, two paired seeds, exactly 400 real unique treatments per target, 400 predictions and disjoint template/reference control pools. Primary references use all right-hand controls; sampled400 is secondary. All candidates share identical truths, templates, references and seeds. Every response fit excludes the entire scored source family.

The initial metadata inventory incorrectly reported two never-count-tested HCT116 targets with 400 cells; the complete older-history check finds only MYRF. MYRF and CAT are reserved; the historical audit below supersedes the initial freshness claim for CAT. H1/HCT development reuses earlier count biology. HEK development prioritizes available fresh targets but is not claimed wholly fresh. HEK confirmation reserves DCAF16 and PSMB9. These four reserved names are disjoint from development, but older V7 split files show CAT, DCAF16 and PSMB9 were previously count-tested. Only MYRF is globally unused. This is reserved algorithm verification on partly reused biological contrasts, not independent blind confirmation. Source summaries were already visible.

A summary screen chooses one of the two interaction penalties. Actual count development compares that model, V13, the fixed shuffled-context control and 1.1x amplification. The plan declares raw-metric guards, a 1% relative expression tolerance and an unanchored six-metric utility. The tolerance is an exploratory practical bound, not a statistically justified noninferiority margin. Promotion requires complete primary panels, improvement in at least four metrics including a DE yield/overlap metric, better utility than V13 and the shuffled control, and aggregate/worst-family guards. Freeze before one confirmation; do not retune after failure. If development fails, retain V13 and leave confirmation unopened.

Control splits balance each source batch between left and right, but family-level truth panels pool treatment batches. This matches the desired common-context origin while retaining potential residual batch-composition confounding. Seeds are repeated draws rather than independent biological lineages. Confirmation covers two families, with no independent H1 confirmation.

## Reproduction

Run from Command Prompt in the project root, with the existing verified source caches and scientific runtime:

```cmd
.\.competition_env\Scripts\python.exe -m unittest test_v14_context -v
.\.competition_env\Scripts\python.exe prepare_v14_experiment.py
.\.competition_env\Scripts\python.exe audit_v14_history.py
.\.competition_env\Scripts\python.exe run_v14_experiment.py screen
.\.competition_env\Scripts\python.exe run_v14_experiment.py oracle --workers 3
.\.competition_env\Scripts\python.exe run_v14_experiment.py development --workers 3
.\.competition_env\Scripts\python.exe assess_v14_experiment.py development
```

Plan creation intentionally refuses to overwrite an existing freeze. Evaluation resumes only hash-matching model and per-case caches. Run confirmation and its assessment only if the frozen development selection promotes a model:

```cmd
.\.competition_env\Scripts\python.exe run_v14_experiment.py confirmation --workers 3
.\.competition_env\Scripts\python.exe assess_v14_experiment.py confirmation
```

Large regenerated model/count-profile caches are ignored by Git. Scripts are in the root; small plans, model sidecars, scores and reports are versioned. Generated counts must be nonnegative integer values and preserve their exact supplied library totals. Tests cover a known context reversal, invalid-value masking, insufficient coverage and held-family label poisoning. No original model, package or frozen evidence is modified. Measured outcomes and the rejection decision are recorded in `V14_RESULTS.md`; no upload is performed by these scripts.

## Metadata correction before selection

The separate `metadata_correction.json` supersedes freshness annotations in the immutable initial plan and split audit. No target/control identities, source masks, seeds, candidates or numerical methods changed. Partial K562/CD4 screen fits were stopped and retained; no numerical screen scores were inspected for selection before the audit. Restarting reuses only hash-verified unchanged fits. Old V7 split filenames were missed by the initial inventory. Consequently reserved verification cannot be presented as independent biological confirmation, regardless of its outcome.

## Completed outcome and final audit

The 72-population development comparison rejected context2. No reserved algorithm verification or submission packaging followed. The nine treatment-oracle cases completed separately and diagnose both transferred-mean and synthetic-distribution limitations. See [V14_RESULTS.md](V14_RESULTS.md) for exact measured values and power limitations.

After evaluation and the once-only development assessment, reproduce the independent split/hash/coordinate/aggregation checks and documentation with:

```cmd
.\.competition_env\Scripts\python.exe audit_v14_experiment.py
.\.competition_env\Scripts\python.exe document_v14_experiment.py
```

The report generator describes this completed rejected experiment and deliberately refuses a promoted selection. `dispatch_v14_oracle.py` was a scheduling-only wrapper while held-family model caches were being built: it waits for completed cache sidecars before calling the unchanged evaluator. Direct sequential execution of the oracle command above produces the same numerical experiment. Its wrapper/evaluator hashes are retained in `oracle_dispatch.json`.
