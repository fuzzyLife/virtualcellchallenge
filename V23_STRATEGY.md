# V23: mean-preserving learned response variation

V22's bounded state-dependent response model improved all three development DE metrics and shape distance, but worsened NMAE beyond its frozen guard. V23 is a separately frozen experiment correcting emitted-mean drift. It retains every V22 factor array and every V16 population response estimate. No historical rejection threshold was changed.

## Integer population correction

Generate the learned conditional population using V22. Generate the exact V16 comparator using the same template cells, libraries and seed. First move counts within rows to match V16's exact gene totals. Then swap counts in 2x2 cell/gene cycles: every cycle preserves both cell libraries and gene totals, while reducing weighted mean-CPM error relative to V16. The code permits at most20,000 cycles and at most the smaller of200,000 or2% of library molecules in cycle transfers. CPM means are approximate, not exactly preserved; residual errors are logged.

This changes per-gene population distributions while protecting averages. No real or hidden treatment values enter the correction. It does not repair a biologically incorrect V16 mean response. A real generated HDAC8 diagnostic reduced mean-probability L1 drift from0.00201614 to0.00000206694 with exact libraries/gene totals; no treatment DE score was computed during that diagnostic.

Frozen development arms: exact V16, original rejected V22 conditional, mean-preserved conditional, equally repaired scalar, equally repaired state-independent, and equally repaired shuffled-state populations. Only mean-preserved conditional is a release finalist. V22 factors exclude the entire evaluated family, including all CD4 conditions when evaluating donor2. Panel genes are excluded from state/program fitting. Basal-state and response-strength proxies remain statistically confounded; this is not identification of individual causal counterfactuals or guide efficiency.

## Validation

Development: the same14 previously studied target/family pairs at two seeds,400 real and predicted cells, six arms =168 comparisons. Existing raw guard thresholds remain fixed; require three improved metrics including DE, improved centered response-projection distance, and superiority to all three equally repaired negative controls on at least two DE metrics. Missing NMAE remains missing. Null predictions must equal V16 bitwise. HCT116's small expression denominator is reported separately; paired absolute error numerators drive its guard.

The Python-heavy swap step scaled poorly across threads. A scheduling-only wrapper runs six isolated family/seed processes, with one case worker and two BLAS threads each. The frozen numeric runner, model arrays, cells, seeds, candidates and guards remain unchanged. Original threaded execution produced no completed cases; its untouched runner remains reproducible. The scheduler collects identical caches and complete panel coverage into the normal analysis files. Source/code hashes and this scheduling amendment are recorded.

## New raw-cell reserve

Metadata-only power checks selected Marson2025 donor2 resting CD4:15 competition targets have >=400 cells, compared with donor1 resting2, donor1 stimulated48h2 and donor4 resting1. The author [Zenodo21729034](https://doi.org/10.5281/zenodo.21729034) file is2,287,711,561 bytes, MD5474833a3f246a748b8076cdcbdcf63e1. Downloaded ranges and assembled file are checked; original donor1 extraction is retained.

Frozen reserve: DEK, DNTTIP1, FXYD5, MED13, MED15, MSI2, PCGF1, PPIP5K2;400 distinct treatment cells each, metadata-selected before reserve scoring. Left/right NTCs are disjoint within23 sequencing lanes (5,006/4,994 controls); all retained cells have low_quality=False. The count matrix is aligned to the competition axis with17,772 measured genes. It represents an existing public study already supplying summary evidence, not a new independent biological study. The reserve model excludes all CD4 raw and summary evidence. Only a locked development finalist opens it once, with the same guard policy. Source-power, lane-mapping and frozen identity amendments are explicit.

## Reproduction

Use the saved plans. Exact scientific runtime: Python3.12.10, NumPy2.5.4, pandas3.0.6, SciPy1.18.1; remaining pins are in requirements-v23-science.txt and analysis/runtime.json. Scientific fitting/scoring/generation/replay use this runtime; the competition environment is used only for unchanged official prep, package validation and lossless decompression. No existing environment dependencies were modified.

Windows CMD:

```bat
cd /d L:\promec\Animesh\Download\Virtual
set "VCC_SCI_PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
"%VCC_SCI_PY%" test_v22_state_response.py
"%VCC_SCI_PY%" test_v23_mean_preservation.py
```

For a new independent experiment directory, prepare V22, fit/run its oracle and model as documented in V22_STRATEGY.md using this scientific runtime. Preserve all existing plans, fits and receipts. V23 preparation refuses to overwrite its frozen plan.

```bat
"%VCC_SCI_PY%" prepare_v23_experiment.py
"%VCC_SCI_PY%" run_v23_parallel_stage.py development --processes 6
"%VCC_SCI_PY%" run_v23_experiment.py null --workers 2
"%VCC_SCI_PY%" assess_v23_experiment.py development
```

Only a passing finalist opens confirmation. The reserve preparation is metadata-only and may run before development finishes. Root collect/probe/extract scripts preserve exact author checksums and original identities. Large data/models/count caches are ignored in Git; small plans, summaries and audit evidence are versioned. Artifact/result documents record the observed final release decision; no leaderboard improvement is guaranteed.
