# V21 results: active variation programs rejected

Completed 2026-10-13. The frozen V21 joint candidate failed development, so no V21 reserve, production model or upload package is released. The cross-version rule retains the separately passing V20 mean model for its own once-only reserve. V16 remains the best published submission.

V21 compares 168 emitted populations across six arms and 28 paired cases. The vocabulary-compatible basis activates every development target group; output coverage is at least .9858. Generated conditional means replay bitwise and all mean-response fields match V20. Thus the rejection is a result of the tested variation mechanism, rather than the old low-coverage abstention.

| Raw family-balanced metric | V16 | Stack mean | Active joint | Matched shuffled joint | Occupancy reweight | Better |
|---|---:|---:|---:|---:|---:|---|
| PDS | .777778 | .766667 | .766667 | .761111 | .766667 | Higher |
| Fidelity | .606476 | .608369 | .608828 | .608251 | .607673 | Higher |
| Reach | .503430 | .517972 | .479652 | .506776 | .504209 | Higher |
| Jaccard | .022365 | .022241 | .022157 | .022276 | .022151 | Higher |
| Eligible-panel NMAE | .849526 | .856223 | .862778 | .853495 | .846275 | Lower |
| Centered program Wasserstein | 2.281830 | 2.284705 | 2.284896 | 2.282672 | 2.286665 | Lower |

The active joint model loses .02378 reach and adds .01325 NMAE versus V16, beyond the .02/.01 guards. Worst-family reach fails too. It improves only two assessed metrics, fails the matched shuffled-program comparator, fails the mean-only DE comparison and worsens centered shape distance versus mean-only. Occupancy reweighting is a diagnostic, not a retrospectively selected release arm.

Independent audits passed: genuine disjoint control/treatment identities, held-family source exclusion, exact prior V16 baseline count hashes, all cache hashes, dense query-state reconstruction, vocabulary and panel exclusion, orthogonal basis, parent mean fields and active norm-matched shuffled programs. Evidence is in `experiments/v21_stack_programs/analysis/development_*`, `selection.json`, and `development_compatibility_independent_audit.json`.

HCT116 NMAE is missing in both seed panels; it was not imputed. The public biology is partly reused, pretraining overlap is unaudited, and these raw metrics do not predict a hidden anchor-scaled leaderboard gain. See V21_STRATEGY.md for the frozen design and commands.
