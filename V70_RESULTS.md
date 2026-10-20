# V70 results: no release

The fixed development screen selected **no eligible candidate**. Raw25 profile-error ratio versus V66 was **1.000641014** (+0.0641%); raw50 was **1.000670064** (+0.0670%). Both lost in Orion; the three development targets in H1 and both CD4 contexts were unchanged. Shuffled raw25/raw50 were worse (1.002234095 / 1.002309382), which did not make the genuine arms improvements.

All six raw sources fitted and the whole-study exclusion audit passed. Synthetic depth-only and known-response checks passed. These establish implementation integrity, not biological transfer accuracy. The frozen >=0.5% profile-gain criterion failed. **No V70 count trial, production model or upload file was generated.** The remaining confirmation target cases were not opened for V70.

The guarded raw25 correction changed only 42/43 CD4 targets, 81 H1 targets, and 64 targets in each Orion line (threshold: maximum absolute LFC change >1e-6, across 300 targets). V71 tests whether preservation of the old model's own geometry blocks useful corrections, alongside batch-adjusted common-program removal and an uncertainty-interval projection. This is a new frozen experiment after V70, not a revision of V70's failed result.

Reproduction: `v70_count_teacher.py freeze`, fit six sources, prepare five public contexts, then `development`; see `V70_STRATEGY.md`. Exact results and hashes are in `experiments/v70_count_teacher/analysis/development_mean_results.json` and `selection.json`. The source audit is `boundary_audit.json`. Production adapters were prepared but were not run. V66 remains best published; no leaderboard gain is claimed.
