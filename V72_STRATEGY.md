# V72: batch-offset pooled response evidence

The first two raw-count teacher pilots failed. Their coverage audit identified a concrete implementation gap: demanding ten treated cells in each sample discarded nearly all Orion targets spread over small sample groups. V72 pools all matched treatment groups with >=1 cell, retaining each group's own NTC composition and total library exposure, before estimating the response and its uncertainty. It still requires >=40 total treated cells and >=20 expected/observed molecules per measured gene.

Observed control proportions define depth offsets. Gamma/quasi-Poisson moments propagate control uncertainty and overdispersion of counts around those offsets. The half-count prior enters posterior uncertainty/response estimation, not the control point exposure. This is an approximate count model, not an exact negative-binomial fit. Own-gene effects remain zero in the learned response and suppression is applied once by the unchanged V66 emitter.

The four V71 arms, bounded corrections, source-study exclusions, development criterion, count confirmation and broad >=30-cell diagnostic are reused without loosening thresholds. The small-batch synthetic test contains 160 treatments spread over 80 batches and would be entirely discarded by the earlier rule; it verifies no spurious depth/batch response and recovery of a known intervention.

Before fitting, inspection caught a temporary binding recursion and clarified the control point estimate. Initial code and plan are preserved in `analysis/implementation_history`; the corrected hashes and explanations are in the plan. No fitted source or outcome existed at those corrections.

Development remains previously used public biology. Selection is locked before the remaining confirmation cases. D3 stays unopened. A production artifact requires generated-profile, worst-study, shuffled and broader diagnostic gates; other failed/missing metrics must remain explicit. None of these proxies guarantees the hidden competition score. No upload or push is automatic.

Use the pinned V24 scientific Python. Run `test_v72_batch_offset.py`, `run_v72_teacher.py` (six source fits, prediction preparation, development), then `v72_batch_offset.py confirmation` only if selected. Run the source/prediction audit, paired-count and broader diagnostic scripts before considering production. Root scripts and compact experiment receipts document all stages.
