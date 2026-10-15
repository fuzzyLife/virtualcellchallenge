# V30: full-target STRING-bilinear and construct-aware response tests

This separately frozen pilot tests two recommendations from POST_V29_RESEARCH.md: broader target-to-response mean prediction and construct-aware source-effect estimation. It does not change V24 or earlier plans/models. Conditional population transport and full pretrained transition adaptation are downstream stages requiring a useful mean finalist; they are not claimed as trained here.

## Data and units

Original Replogle2022 K562 GWPS and RPE1 pseudobulk files are already local and publisher-checksummed. Means are fractional per-cell count summaries, not raw integer single cells. Use all qualifying targets (>=30 filtered cells; observed source fold_expr<=0.5), aggregate duplicate output symbols, and normalize on the common competition-measured source axis BEFORE filtering baseline expression. This common-axis normalization is a source-transfer convention; it is not identical to full-transcriptome source normalization or the competition count contract.

The source inventories contain8,443 K562 and1,891 RPE1 targets. STRING coverage leaves8,110 and1,889 fitting targets, respectively. RPE1's zero direct panel overlap does not prevent auxiliary structured learning. H1/Orion/CD4 summaries use the existing validated response cache and its measured masks. Those cached target rows receive equal weights; their constant n=100 storage value is only a weight placeholder, cancels at within-source weight normalization, and must not be interpreted as a measured cell count. Cell-count-based weights are used for the original Replogle pseudobulk targets.

Use the author's pinned SLIM STRING HDF5 embeddings, first64 coordinates as in its default loader. Revision5a7e9ade5d0a6b6331e6dbc81181450605047bcc; Git blob8282fd198c2713f6d0a823c88211b6d76308a7df; SHA256789416877b8701ef6f800106d26bf7bb97ea8e72744e6ab93e24933a717f247d. The artifact's19,502,608 bytes and Git blob were verified. Author software MIT license and receipt are retained beside the local data. Cite Hu et al., SLIM, https://doi.org/10.64898/2026.08.07.743481 and its STRING embedding reference https://doi.org/10.1038/s41540-026-00746-8.

## Model

Fit a separate training-only response basis for each source, then weighted bilinear ridge from64 STRING coordinates into16/32 source response programs. This is an independently implemented SLIM-inspired source-transfer adaptation, not a reproduction of its benchmark. Target features are centered/scaled using that source's fitting targets only. The bilinear loss excludes each perturbation's own output exactly. Own coordinates receive centered-zero initialization for SVD only; conjugate-gradient normal equations remove their model terms and labels. An independent explicit masked least-squares fixture checks the solver. Save bases, weights, scales, predictions and original fitting target names for replay.

Pool only available source outputs. Query-control rank compatibility weights sources; Replogle K562/RPE1 share one family vote and CD4 conditions share another. Available-source counts are normalized separately per target/output coordinate, so a missing measurement abstains. Missing target embeddings also abstain. Correct V24 by alpha0.125/0.25 times a masked difference bounded at+/-0.5log2FC. This deliberately tests response means with the prior emitter held for any later count phase.

Variants: full-target STRING, panel-only STRING, randomized embeddings, and construct-aware full-target STRING. Ranks16/32 and lambda1/10 are screened; ridge penalty is lambda*n_fitting_targets/64. Where fewer than20 source targets are eligible, explicitly abstain rather than fitting or inventing labels. In particular panel-only RPE1 abstains.

## Construct-aware teacher and its limits

For targets with at least two effective published constructs, correct each construct response using observed source suppression with gain clipped0.8-1.8. Estimate between-construct heterogeneity after a bounded Poisson uncertainty proxy; pool with precision weights capped at10 times the minimum. Preserve heterogeneity rather than dividing it away as replication noise. Single-construct targets retain ordinary aggregation. There are172 K562 and37 RPE1 targets with multiple effective constructs.

Promoter/transcript constructs are distinct interventions, not technical replicates or individually resolved sgRNA efficiencies. This pilot does not identify guide efficiency from anonymous query controls and does not fit unsupported guide random effects in H1/Orion composite-guide labels. It is the feasible construct-aware part of the broader guide/assay recommendation; author RNA-fingerprinting and a full guide-resolved hierarchical teacher remain separate future work.

## Prespecified evaluation

Nested mean selection treats H1, Orion (HCT116 and HEK293T together), and CD4 (all conditions together) as outer families. For each outer family, select rank/lambda/blend using only other inner families while excluding BOTH outer and inner from every predicted source. Source-specific bases are never mixed across families; source models trained for other evaluations may exist on disk but excluded source models are never requested in that query.

The primary mean screen uses all available masked summary-coordinate LFC errors, excluding own genes. It is an exploratory screen on reused public biology, not blind validation or a private anchored competition score. Equal family weighting prevents CD4 condition count dominating. A primary candidate must improve MSE>=1% and MAE>=0.5%, stay within2% of baseline MSE in its worst family, and beat matched randomized/panel-only controls on both error measures. Construct controls use the construct candidate's chosen hyperparameters. Selection must be recomputed independently.

Only a passing mean candidate proceeds to prespecified400-cell count/null tests. Only a passing count finalist can open the still-untouched donor reserve, with all CD4 evidence excluded during its evaluation. There is no automatic package or upload.

## Freeze record and commands

Final plan SHA256:5c25fb3c10ded81c91f47900c29500957dbcd322cbcdaf62693ae38984c8c549. The initial plan is preserved; a documented pre-outcome amendment adds explicit zero-overlap abstention and coordinate-wise available-source family normalization. It occurred before any model was trained or any candidate outcome opened. Earlier scientific code/results are unchanged.

Use the scientific runtime recorded in the plan (Python3.12.10), not the separate official packaging runtime:

```bat
"C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe" test_v30_bilinear.py
"C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe" run_v30_experiment.py prepare
"C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe" run_v30_experiment.py train
"C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe" run_v30_experiment.py summary
"C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe" audit_v30_experiment.py
```

`freeze` is used once for a new untouched experiment directory; do not overwrite an existing plan. Four worker processes use two BLAS threads each. Large matrices, embeddings and models remain outside Git; source code, hashes, source/model metadata, plans and small evidence are versioned. The independent auditor recomputes source profiles and selected original target summaries, checks fitted masked-loss gradients, replays predictions and verifies actual excluded-family model requests.
