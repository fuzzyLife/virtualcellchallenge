# V58: held-study guard for the cross-dataset PIE experiment

This corrects an implementation dependency found while reviewing V57, before count generation. The old V57 plan, inputs, neural predictions, calibrators and any response results are preserved. V57 is not evidence of clean CD4 study exclusion: although its regression labels used only H1/Orion, baseline features for those source examples could include held-CD4 response-bank values.

The corrected baseline excludes **CD4 plus the source recipient study** in every calibration example. H1 excludes H1+CD4. HCT116/HEK293T exclude both Orion lines+CD4. CD4 excludes all CD4. An independent metamorphic check replaces every CD4 bank response and error with large constants and requires every baseline-bank field to remain bitwise unchanged. Matching query responses are still removed from PIE donor evidence as before.

The neural forward passes remain valid: baseline is not an input to PIE. Preparation checks all neural input fields, all response labels, gene axes, masks and query controls bitwise against V57, then copies the unchanged three-head outputs with their parent SHA256 receipts. The CD4 comparator baseline must also remain bitwise unchanged.

Numerical calibration/decoder bodies are identical to frozen V57 after rebinding the experiment imports. Architecture, features, ridge grid, whole-target split, samples, count cohorts, generation seeds and 176-count-comparison design remain fixed. This is a leakage repair, not a choice based on CD4 accuracy. The public biology was already used in earlier project development; D3 stays unopened.

The original PIE evidence and assets remain under the V57 cache and audited receipts; V58 stores only corrected inputs, calibrators and results in its own namespace. One masked NCBI alias remains explicitly documented for TMEM104; no mixture of neural configurations or post-result feature tuning is introduced.

Commands, using the existing Python 3.12 scientific runtime:

```bat
set PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe
%PY% -X utf8 prepare_v58_pie.py freeze
%PY% -X utf8 run_v58_response.py
rem Inspect response_summary.csv before generating counts.
%PY% -X utf8 run_v58_count_pipeline.py
%PY% -X utf8 export_v58_pie_arrays.py
```

The freeze refuses to overwrite an existing plan. Root scientific scripts and compact receipts are versioned; regenerable arrays are excluded from Git. Read-only evidence caching is a performance optimization checked for exact numerical parity. No scoring API upload occurs. Model and code details, the evidence projection and dataset overlap limitations are documented in the preserved [V57 protocol](V57_PIE_CROSS_STRATEGY.md) and [V56 audit](V56_PIE_RELEASE_REVIEW.md).
