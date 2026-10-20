# V78 results: exact means protected; shape restoration did not qualify

All110 development populations passed the independent raw-statistic audit: every baseline was regenerated from V66, raw profiles were bitwise identical, all BH arrays and four DE metrics replayed, PDS recomputed, and880 tied-rank coordinates checked. This establishes integrity and isolates cell distribution from bulk response quality.

| Arm | Raw fidelity | Raw Jaccard | Raw reach |
|---|---:|---:|---:|
| V66 |0.579503341|0.062351815|0.214770418|
| control10 |0.579138900|0.062321119|0.206051567|
| control25 |0.580175598|0.062687425|0.213708847|
| learned25 |0.579584994|0.062035321|0.215831990|
| query50 |0.484711488|0.040880156|0.226967777|
| query100 |0.384827814|0.014166104|0.289934209|

PDS is exactly0.833333333 for all arms on this small development panel. No sample returned the fallback. The available three-study NMAE means remain close to baseline0.932783543. These are unscaled public metrics and cannot be compared numerically to hosted anchor-scaled component scores.

Control25 has the most favorable small precision/overlap change, but it misses both prospective minimum gains: fidelity+0.000672257 versus required+0.002, Jaccard+0.000335610 versus required+0.0005. It was not promoted by lowering those cutoffs. Larger restoration materially harms fidelity and overlap despite raising reach; the different DE objectives conflict.

No arm was eligible. No remaining-target confirmation, repeated seed or V78 production package was generated. Full control variance is not established to equal treated-cell variance, and projecting proposals onto fixed means changes their marginal distributions. These failures reject the executed proposal; they do not prove that cell heterogeneity is irrelevant.

See `MARGIN_VARIANCE_STRATEGY.md`, `v78_margin_variance.py`, `run_v78_variance.py`, and the complete result/audit files in `experiments/v78_margin_variance/analysis`. The exact original auditor source is retained with its recorded hash in `auditor_history` after adding support for later versions.
