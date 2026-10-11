# V13 published leaderboard result

Published as `AvAsYkuxFSjbxYDCd0wJ`, overall **0.137883134**, a gain of **+0.000879160** (0.64%) over V11 on the same panel and anchors. This is the best recorded project validation score. Rank was **558** in the user report and **559** at the API check; ranks change with other entries.

## Scaled comparison

| metric | v11 | v13 | change_from_v11 | overall_contribution_change_from_v11 |
| --- | --- | --- | --- | --- |
| overall | 0.137004 | 0.137883 | 0.000879 | undefined |
| pds | 0.638267 | 0.639136 | 0.000869 | 0.000145 |
| mse | 0.024946 | 0.032878 | 0.007932 | 0.001322 |
| nmae | 0.081058 | 0.081822 | 0.000763 | 0.000127 |
| fid | -0.015859 | -0.015618 | 0.000241 | 0.000040 |
| reach | 0.080022 | 0.075574 | -0.004448 | -0.000741 |
| jac | 0.013590 | 0.013508 | -0.000082 | -0.000014 |

All displayed components are scaled scores: higher is better, including MSE and NMAE. Four components improve. Expression accuracy contributes +0.001322 to overall, partly offset by reach (-0.000741). Contributions are component changes divided by six; they reconcile exactly with the overall gain and do not establish causal attribution.

## Raw comparison

| metric | v11 | v13 | change_from_v11 | relative_change | improved_from_v11 |
| --- | --- | --- | --- | --- | --- |
| pds_cosine | 0.789049 | 0.789487 | 0.000438 | 0.000556 | True |
| expr_mse_unbiased_capped_norm | 0.987167 | 0.976550 | -0.010617 | -0.010755 | True |
| de_wilcoxon_lfc_nmae | 0.951606 | 0.951212 | -0.000393 | -0.000414 | True |
| de_wilcoxon_direction_fidelity_yield_raw | 0.508011 | 0.508106 | 0.000095 | 0.000187 | True |
| de_wilcoxon_direction_reach_raw | 0.149502 | 0.145631 | -0.003871 | -0.025893 | False |
| de_wilcoxon_sig_jaccard | 0.035282 | 0.035250 | -0.000032 | -0.000914 | False |

Raw expression error decreases **1.08%**; raw reach decreases **2.59%**. Fidelity is still weak and its scaled component remains negative. These results support a small calibration improvement, not a broad solution to context transfer.

## Interpretation and next tests

V13 corrects unquantifiable on-target RNA being amplified as weak measured knockdown. It also includes small pooling-precision differences recorded before publication. Its population emitter and seed match V11. A single published result does not isolate these effects, establish statistical significance, or demonstrate performance on unseen final-round contexts. The frozen local count-confirmation result remains inconclusive; do not rewrite it as a pass.

Preserve V13 as the best recorded validation-score reference and V11 as a comparison with better reach and overlap. Prioritize a separately frozen test of the published construct-level K562 knockdown metadata, and coverage guards that preserve V11's reach. Neutral-population false-DE bias and unrealistic dispersion remain unresolved. Use held-family validation before another upload; this score alone does not justify tuning a larger model.

## Provenance and reproduction

Submitted `2026-10-11T07:28:25.755195Z`. Upload took 16 minutes 27 seconds in the user report. Partition `val`, panel `vcc2026-val-1`, anchors `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`. Package SHA-256: `2454a8b8ce694217962b89f1bac4ec0fb1bbd01a396c1ffe5cd780f4baa11591`. The status API does not expose the uploaded hash; the user's upload record links the local audited artifact to this entry.

All frozen code, selection, confirmation, generation and integrity records remain unchanged, including generation-time `uploaded=false`. The original report is archived byte-for-byte. The release manifest still describes the generation-time files; `publication_manifest.json` maps its original `V13_SUBMISSION.md` hash to the archive and separately hashes current publication documents. Rerun this recorder after the original document generator to restore the publication overlay.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc status AvAsYkuxFSjbxYDCd0wJ
.\.competition_env\Scripts\python.exe record_v13_result.py --api-file experiments/v13_response_programs/analysis/published_api_status.json
```
