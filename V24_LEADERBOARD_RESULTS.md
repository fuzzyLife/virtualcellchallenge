# V24 published leaderboard result

V24 is published as `81Y2yX7utdAau9rr3jCj`, overall **0.142997855**, improving on V16 by **+0.001383998** (0.98%) and V20 by **+0.003986288**. It is the new best recorded project validation submission. All six scaled components improve versus V16.

Rank was **530** at the user report and **531** at the read-only API check. Historical ranks change as other teams submit; the same-panel score comparison establishes the gain.

## Scaled scores

| metric | v16 | v20 | v24 | change_from_v16 | change_from_v20 | overall_contribution_vs_v16 |
| --- | --- | --- | --- | --- | --- | --- |
| overall | 0.141613857 | 0.139011567 | 0.142997855 | 0.001383998 | 0.003986288 |  |
| pds | 0.641552935 | 0.639288077 | 0.644246367 | 0.002693432 | 0.004958289 | 0.000448905 |
| mse | 0.046085072 | 0.047823003 | 0.047568407 | 0.001483336 | -0.000254596 | 0.000247223 |
| nmae | 0.080344812 | 0.079577327 | 0.082270573 | 0.001925760 | 0.002693246 | 0.000320960 |
| fid | -0.017295588 | -0.016123050 | -0.015791837 | 0.001503750 | 0.000331213 | 0.000250625 |
| reach | 0.084381077 | 0.068935293 | 0.084839840 | 0.000458763 | 0.015904547 | 0.000076460 |
| jac | 0.014614833 | 0.014568752 | 0.014853782 | 0.000238949 | 0.000285030 | 0.000039825 |

All six displayed component scores are higher-is-better, including expression MSE and LFC NMAE. Overall is their unweighted average; the component contributions above sum exactly to the overall gain. This is arithmetic decomposition, not causal attribution. The improvement is consistent across components but modest, not a ranking breakthrough.

## Raw metrics

| metric | v16 | v24 | change_from_v16 | relative_change | higher_is_better | improved |
| --- | --- | --- | --- | --- | --- | --- |
| pds_cosine | 0.790561130 | 0.791865478 | 0.001304348 | 0.001649901 | True | True |
| expr_mse_unbiased_capped_norm | 0.961835357 | 0.960414881 | -0.001420477 | -0.001476840 | False | True |
| de_wilcoxon_lfc_nmae | 0.952086028 | 0.950685826 | -0.001400203 | -0.001470668 | False | True |
| de_wilcoxon_direction_fidelity_yield_raw | 0.507558266 | 0.508049590 | 0.000491323 | 0.000968014 | True | True |
| de_wilcoxon_direction_reach_raw | 0.153714413 | 0.154141828 | 0.000427415 | 0.002780577 | True | True |
| de_wilcoxon_sig_jaccard | 0.035660215 | 0.035753233 | 0.000093018 | 0.002608454 | True | True |

Raw expression MSE and raw LFC NMAE are lower-is-better; other raw metrics are higher-is-better. These values are not the CLI's scaled scores. Direction fidelity remains below its zero-score anchor; low reach and significance overlap remain substantial gaps.

## What this supports

The final expert changes 207 context-A target groups, with maximum coordinate correction0.25 log2FC, compatibility and donor-direction confidence gates, and condition weights from measured controls. The remaining693 groups reuse exactly verified V16 count blocks; original control sampling, libraries and the residual10 emitter are retained. This supports a narrow source/context calibration over another broad generic ensemble. It does not isolate the effect of individual genes, condition weights or donor gates.

Local reserve fidelity, discrimination, expression and LFC accuracy improved, while reach and overlap dipped slightly. All six hosted components now improve versus V16. The host outcome does not retroactively remove public-study overlap, missing development NMAE, shared upstream nuisance preprocessing or limited reserve-target coverage. V22 and V23 remain rejected experiments. Larger causal population-variation improvements remain unresolved.

## Provenance

Partition `val`, panel `vcc2026-val-1`, anchors `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`; same as V16/V20. Submitted `2026-10-14T17:25:19.507408Z`. Upload took16 minutes18 seconds. Package is4,073,656,320 bytes; SHA256 `42cd3906e17aa2e8e546efea44f1e726a38e577605f211ad96f1564362a20cc6`. Official prep, all2,368,814,375 packaged count/index comparisons, every control identity/library check and12 exact group replays passed before submission.

The API confirms the entry and scores but exposes no uploaded-file hash. The user's upload command connects the audited local package to this entry. Generation-time uploaded=false snapshots, numerical code, model arrays and packages remain unchanged. Original report documents are archived byte-for-byte; publication_manifest.json maps generation seals and current publication evidence. Re-running the report generator restores this publication overlay from saved API evidence.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc status 81Y2yX7utdAau9rr3jCj
.\.competition_env\Scripts\python.exe record_v24_result.py --api-file experiments/v24_cd4_response/analysis/published_api_status.json
```
