# After V82: why the plateau needs an evidence test

The best published project score remains V66/V82, **0.1508821875498449**. The response-transfer prototype and its completed decision are documented in `V84_STRATEGY.md` and `V84_RESULTS.md`. This review distinguishes observed project limitations from hypotheses. It does not establish what private leading teams are doing.

## What the project evidence says

1. **Much of the search has reused the same response information.** The current response source cache still comes from K562, CD4, H1, HCT116 and HEK293T families. V10's population transport was a large improvement; V63 identity protection and V66 strength calibration were smaller gains. Many other versions alter mixtures, means or synthetic count realization around those sources. Different version numbers do not imply new experimental evidence.
2. **Local improvements have repeatedly failed to transfer.** V69 regressed from 0.150882 to 0.095263 on the same hosted panel. V70-V81 then tested twelve families with 410 count comparisons; none met its original public requirements. These extensive checks reuse public biology and cannot provide 410 independent confirmations. V82 reproduced V66 exactly, confirming that the recovery score itself was stable.
3. **New context labels can hide zero target coverage.** The pinned Arc essential-screen DE release has no direct target overlap with our 300 genes in K562, RPE1, HepG2 or Jurkat. The empty restricted HepG2/Jurkat tables are correct abstentions, not an ingestion defect. V17 already recorded this restriction; repeating direct augmentation without first carrying that fact forward was a planning mistake. The reason for publisher target absence is not established, and it does not imply that all competition genes are nonessential.
4. **Cell similarity alone does not establish response similarity.** Basal rank correlations can route available evidence, but they do not reveal how an unseen perturbation changes regulation. This is a limitation of our transfer assumptions, not proof that no basal information can help.
5. **The submitted object is a count population.** Intended logFC, realized bulk profiles, significance and directional calls can move differently. V69 and the later fixed-margin experiments show why tuning a response array or matching its mean is insufficient evidence for improving all six hosted metrics.

These support an inference: the next substantial gain is more likely to require transferable intervention information, or a materially better way to identify it, than another blend of the same submitted responses. They do not prove a universal accuracy ceiling.

## Relevant publications and what can actually transfer

| Primary source | Relevant idea | Limitation for this competition |
| --- | --- | --- |
| [Systema, Nature Biotechnology](https://www.nature.com/articles/s41587-025-02777-8) | Common perturbation-related variation can make models look accurate without reconstructing intervention identity. Include strong constant-response comparators. | Its evaluation framework is not the competition's six-metric scoring system; use it as a diagnostic rather than substitute its score. |
| [scDEBART, ICML 2026](https://proceedings.mlr.press/v306/sung26a.html) | Learn basal-conditioned expression changes rather than only reconstruct static expression. | The reported results include fine-tuning on Perturb-seq datasets. They do not establish zero-shot performance in our controls-only target contexts. |
| [C3TL, primary methods](https://arxiv.org/html/2603.13051v1) | Share representations across measured intervention/context pairs. | Its context encoder aggregates measured effects in that context. A missing intervention/context pair is a different information setting from having no target-context intervention responses. |
| [Virtual Cells Need Context, Not Just Scale, preprint record](https://pmc.ncbi.nlm.nih.gov/articles/PMC12919078/) | Context diversity and transport are important questions to test, rather than assuming scale resolves them. | A preprint and a motivation for experiments, not an independently verified VCC recipe. Subsequent automated reopening of the record was access-limited; no new full-text claims are inferred from that failure. |

The scDEBART abstract, author preprocessing/training code and release inventories were inspected; its full PDF could not be retrieved in this run. No pretrained checkpoint was downloaded or executed. The C3TL HTML methods were inspected. The response-transfer model implemented here is our own reduced-rank prototype, not a reproduction of either paper or their reported results.

## Released code also needs verification

The inspected scDEBART author branch is pinned to commit `8996c7ca8c336d636de4ba06a4d5d97b2bbe2ce8`; its dataset inventory is pinned to `3d9ba7f3eeb639675f5e6bf44ce61c8cc707dd5a`. Apache 2.0 license text and the reviewed code are retained with their hashes.

`audit_v84_paper_code.py` extracts only the two reviewed attention class definitions and tests their mask behavior under our PyTorch 2.10.0 environment. All-false boolean masks in these classes block encoder/cross attention; the inspected triangular boolean decoder mask permits later-token influence. The inspected training parser also lacks several architecture arguments used by the model. Numeric fixtures and exact source hashes are in `scdebart_code_audit.json`.

This establishes behavior of that inspected code in this environment. It does **not** establish which code trained the released checkpoints, invalidate the paper's results, or show that the method cannot work. A faithful checkpoint evaluation would require resolving its architecture and training provenance before quoting performance. Silently patching masks and calling it the published model would be misleading.

## The experiment performed instead

V83 first tried direct new-context response augmentation and stopped at the zero-coverage check. V84 then asks a distinct question: can **other shared interventions teach how a measured K562 response changes with context**?

Four SHA-verified source tables supply 638 paired interventions on 6,087 shared nonpanel outputs. A fixed 513/125 intervention split trains a rank-32 response-difference operator and a simple diagonal comparator. We compare genuine fitting pairs against refitted deranged pairs, copying K562, and a constant destination response. Transfer corrections to our existing target responses use only NTC expression for context routing.

Evaluate both the 655 reused public pairs and a newly extracted D3 donor from a previously used study. The donor contributes 14 metadata-qualified targets, disjoint fit/reference NTCs, matched 400-cell populations, and two technical seeds. Every CD4 response source is excluded. No donor treatment outcome is used for fitting or choosing the model.

If the primary fails, a saved conditional diagnostic supplies treatment information **only to non-deployable oracle populations**. It asks whether a more realistic population shape helps at the same mean, versus whether the measured fitting-half mean is needed. This can identify the next bottleneck without turning donor labels into apparent zero-shot accuracy.

No additional upload is justified solely by novelty, source-task accuracy, or a favorable partial metric. The actual development and donor outcomes, failures and scope are recorded together in `V84_RESULTS.md`.

## Diagnostic conclusion after the tests

The fixed response operator failed public development and the first donor seed. The presaved six-target oracle then supplied measured fitting-half means to the existing population model. Profile error fell by 34.72%, and raw directional fidelity rose from 0.4981 to 0.8287. This establishes headroom in those cases; it does not show that our available prediction inputs can recover the missing means.

Using real fitting-half cells instead of the old shape increased DE Jaccard from 0.1608 to 0.3145 at identical bulk profiles. Distribution also matters. The projected fitted-shape-at-baseline-means branch failed its CPM contract, so its poor scores do not constitute a clean shape-only result. Disjoint fitting/scoring halves prevent a trivial same-cell reconstruction, but both halves belong to the same donor and the sample is only six targets.

The priority is to demonstrate target-specific response-mean improvements in the actual scoring statistics, on a genuinely separate study, before another count generator. New measured evidence must pass panel coverage checks or a genuine paired-intervention transfer test that beats both a constant response and a refitted shuffled control. A new architecture name or another ensemble is insufficient evidence. V84 tested this route and did not meet its original release criteria; the complete second-seed outcome is included in the results report.
