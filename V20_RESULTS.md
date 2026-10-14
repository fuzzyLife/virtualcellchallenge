# V20 development results

Completed development on 2026-10-13. The mean ensemble passes the frozen available-metric guards and negative controls. Its development gains are small and mixed; a reserved-case result is required before release. The joint and reweight arms in this original experiment abstain on every case because the inherited program basis has insufficient Stack vocabulary coverage, so their emitted counts are exactly the mean arm. V21 separately fixes compatibility and rejects the active joint mechanism.

| Raw family-balanced metric | V16 | Stack mean | Shuffled response | Amplitude control | Uncalibrated mean | Better |
|---|---:|---:|---:|---:|---:|---|
| PDS | .777778 | .766667 | .761111 | .777778 | .597222 | Higher |
| Fidelity | .606476 | .608369 | .605512 | .606630 | .511648 | Higher |
| Reach | .503430 | .517972 | .464301 | .452599 | .153912 | Higher |
| Jaccard | .022365 | .022241 | .022260 | .022140 | .018586 | Higher |
| Eligible-panel NMAE | .849526 | .856223 | .856070 | .856473 | 2.448962 | Lower |

Mean ensemble changes: fidelity +.001894, reach +.014542, PDS -.011111, Jaccard -.000124, NMAE +.006697. Paired relative expression numerator improves only .0178%. Exactly three assessed metrics improve: fidelity, reach and expression. It beats amplitude on three DE measures and shuffled response on two. All predeclared available overall/worst-family tolerances pass. These are not anchor-scaled competition scores. Raw expression ratios have small reference denominators; the selection uses the paired relative numerator change, not the displayed ratio alone.

Independent raw-input reconstruction and case/cache audits pass. All 28 V16 count digests match the previously tested V17 comparator exactly. Checkpoint/code hashes, strict full loading, original-author mean/dispersion bitwise parity and deterministic repeat checks pass. HCT116 NMAE is unavailable in both development seed panels and is not imputed. Public biology is partly reused and pretrained study overlap remains unaudited.

Evidence: `experiments/v20_stack_transfer/analysis/plan.json`, `selection.json`, `development_scores.csv`, `development_panels.csv`, `development_summary.csv`, `development_independent_audit.json`, raw-input audit, `mean_adapter_parity.json` and pinned asset/runtime manifests. See V20_STRATEGY.md for formulas and reproduction. The final reserved result and release decision are appended after the once-only assessment completes.

## Reserved verification: passed

The fixed mean candidate passed the once-only reserve on 2026-10-13: **80 complete-count comparisons**, eight context/target pairs and two seeds, including the fixed V16, amplitude, shuffled and inactive-reweight controls. The final choice was recorded in `reserve_entry_decision.json` before reserved scoring. No coefficients or source rules were changed after seeing this result.

| Raw family-balanced metric | V16 | Stack mean | Change | Better |
|---|---:|---:|---:|---|
| PDS | .805556 | .805556 | .000000 | Higher |
| Fidelity | .606165 | .607219 | +.001054 | Higher |
| Reach | .514811 | .536862 | +.022051 | Higher |
| Jaccard | .028855 | .028754 | -.000101 | Higher |
| Eligible-panel NMAE | .852138 | .847810 | -.004328 | Lower |

Paired relative expression numerator improves **.2718%**. Fidelity, reach, eligible NMAE and expression improve; PDS is unchanged and Jaccard declines slightly. Every available overall/worst-family guard passes. The candidate beats amplitude on two DE measures and shuffled response on three. Centered program Wasserstein decreases from 2.374581 to 2.371434, a secondary observation rather than an active-variation claim.

Independent case/cache and reserved source-input audits pass. HCT116 NMAE is still unavailable in both seed panels. Consequently this is **available-metric exploratory public algorithm verification**, not complete six-metric or independent biological confirmation. Production is allowed by the frozen rules, but a private leaderboard gain remains unknown.

Reserved evidence: `confirmation_result.json`, `confirmation_evaluation.json`, `confirmation_summary.csv`, `confirmation_panels.csv`, `confirmation_scores.csv`, `confirmation_expert_diagnostics.csv`, `confirmation_independent_audit.json`, `confirmation_input_preparation.json`, and `confirmation_source_input_independent_audit.json` under the same analysis directory.

## Production artifact: complete and independently validated

The final A/B/C matrix contains 360,000 cells, 900 groups and 18,533 genes. **864 groups** receive a Stack mean correction; **36** bitwise-unchanged groups reuse the independently verified V16 output. The largest pre-emission correction is **0.075611755 log2**, below the frozen .125 bound. All twelve unsupported targets retain V16 in each context: C5orf22, CDK8, GLYR1, IL10RB, MAN2A1, RFK, SCAF8, SPATA20, STT3A, TMEM147, UNK and UROS.

The 24-group count preflight passed. Official VCC prep passed, and the independent validator reconstructed all 900 correction vectors, compared every one of **2,368,754,215 packaged nonzero values** plus indices/pointers/axes, verified all original cell/guide/library identities and exactly regenerated 12 groups. The separate input audit checked 1,200 genuine query controls and 233,296 source count rows across 587 source-target prompts, including reused conditioning controls. Metadata-only export and runtime compatibility fixes are archived and documented; no selected parameters or reserved outcomes were changed.

The upload package is **4,073,574,400 bytes**, SHA-256 `8370080f290febfbec59e25a54b175317f1a22609ee3c90faa1997774611ac50`. It has not been uploaded by this workflow and has no new leaderboard score. See [V20_SUBMISSION.md](V20_SUBMISSION.md) for the artifact and Windows CMD upload command.

## Published result

Published `B92qai3cUvEHeobzp8Fk`, overall **0.139011567**, below V16 by **0.002602290** (1.84%). Expression and fidelity improve; reach and three other scaled components decline. V16 remains the best model. See [V20_LEADERBOARD_RESULTS.md](V20_LEADERBOARD_RESULTS.md) for exact comparisons and limitations. The preceding generation report is archived unchanged.
