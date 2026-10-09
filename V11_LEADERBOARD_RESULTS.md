# V11 published leaderboard result

Published as `bmIAGveuEvFiYQTvz0Ma`, overall **0.137003974**, rank **529** in the user report and **529** at the recorded API check. V11 is the best published project submission, improving on V10 by **+0.020249853**. Ranks change over time. V4, V9, V10 and V11 use the same validation panel and anchors.

## Scaled component comparison

| metric | v4 | v9 | v10 | v11 | change_from_v10 | overall_contribution_change_from_v10 |
| --- | --- | --- | --- | --- | --- | --- |
| overall | -0.050674 | -0.070447 | 0.116754 | 0.137004 | 0.020250 | undefined |
| pds | 0.255302 | 0.313301 | 0.572464 | 0.638267 | 0.065803 | 0.010967 |
| mse | 0.000000 | 0.033719 | 0.049815 | 0.024946 | -0.024868 | -0.004145 |
| nmae | 0.067539 | 0.074927 | 0.049886 | 0.081058 | 0.031172 | 0.005195 |
| fid | -0.683855 | -0.921183 | -0.022335 | -0.015859 | 0.006475 | 0.001079 |
| reach | 0.086849 | 0.102034 | 0.046547 | 0.080022 | 0.033475 | 0.005579 |
| jac | -0.029878 | -0.025482 | 0.004148 | 0.013590 | 0.009441 | 0.001574 |

Five of six components improve over V10. The expression-accuracy component declines. Scaled MSE and NMAE are accuracy scores: higher is better. Component contributions are changes divided by six and sum to the overall change; they do not isolate the effects of individual model changes.

## Raw metric comparison

| metric | v4 | v10 | v11 | change_from_v10 | improved_from_v10 |
| --- | --- | --- | --- | --- | --- |
| pds_cosine | 0.616462 | 0.759535 | 0.789049 | 0.029513 | True |
| expr_mse_unbiased_capped_norm | 1.247733 | 0.945955 | 0.987167 | 0.041212 | False |
| de_wilcoxon_lfc_nmae | 0.959164 | 0.970311 | 0.951606 | -0.018706 | True |
| de_wilcoxon_direction_fidelity_yield_raw | 0.309235 | 0.505983 | 0.508011 | 0.002028 | True |
| de_wilcoxon_direction_reach_raw | 0.156391 | 0.120318 | 0.149502 | 0.029184 | True |
| de_wilcoxon_sig_jaccard | 0.019731 | 0.031815 | 0.035282 | 0.003467 | True |

Raw reach increases **24.26%**, fold-change error decreases **1.93%**, and expression error increases **4.36%**. Fidelity yield improves slightly but its scaled score remains below the anchor. Reach remains below V4's raw reach, despite the higher overall score.

## What we learned

V11 improves on V10 on the official validation panel. The two changes are consensus-supported recovery of shrunken source responses and increased synthetic-cell dispersion. Their combined result improves discrimination, fold-change accuracy, fidelity, reach and overlap, while worsening expression error. The public confirmation checks anticipated the average reach/fold-change gains and the expression-error tradeoff. A single submission cannot isolate which change caused each improvement or establish performance in unseen final-round contexts.

Preserve V11 as the new published baseline. The next controlled experiments should protect its direction and coverage gains while improving expression accuracy: separately ablate response rescue and dispersion, then test gene-specific uncertainty/dispersion or mean-versus-pooled response calibration. Choose any new parameters using held-family validation and a separate confirmation set, rather than tuning against this score alone. The substantial neutral-control false-DE bias and SLIRP power-check reach loss remain unresolved; a better leaderboard score does not establish realistic biological cell variance.

## Provenance and reproduction

Entry submitted `2026-10-09T08:16:34.765565Z`. User-reported upload duration: 16 minutes 36 seconds. Partition `val`, panel `vcc2026-val-1`, anchors `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`. Local independently audited package SHA-256: `7841429ef987c43abcfb01a5be8f09864c5f8955e1bd7aae68fb4a03240be31d`. The status API does not expose the uploaded file hash; the local artifact and server publication are linked by the user's upload record, not an independently confirmed server-side hash.

All frozen numerical code, generation choices, confirmation and count-audit records remain unchanged. The original submission document is archived byte-for-byte as `experiments/v11_consensus_coverage/analysis/generation_submission_document.md`, matching its frozen documentation-validation hash. `publication_manifest.json` links the publication records and current documents separately.

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc status bmIAGveuEvFiYQTvz0Ma
.\.competition_env\Scripts\python.exe record_v11_result.py --api-file experiments/v11_consensus_coverage/analysis/published_api_status.json
```

These commands inspect or document the existing entry. If the generation pipeline regenerates the pre-publication report, rerun the recorder to restore the publication overlay.
