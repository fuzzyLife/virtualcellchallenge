# V71 results: no release

The fixed development criterion selected no candidate. Equal-study desired-profile error ratios versus V66 were raw25 **1.005377339**, center25 **1.001455043**, center50 **1.006805253**, and credible25 **0.999990482**. The last is only a 0.00095% decrease, far below the prespecified 0.5% minimum. No count trial, production prediction or submission was generated. The confirmation cases remain unopened.

Raw25 without the identity guard was worse than guarded V70 raw25, so this test does not establish that removing the guard is beneficial. The audit then exposed a more concrete limitation: V70 counted targets as eligible even where its per-batch minimum left zero usable evidence. Only 12/25 H1, 13/275 HCT116 and 23/287 HEK293T targets had supported coordinates; CD4 had support for all listed targets. V72 corrects that data-aggregation rule, preserving both failed experiments.

All values are repeated-public-development diagnostics, not leaderboard forecasts. V66 remains highest published. See `development_mean_results.json`, `selection.json`, the original source metadata and the frozen V71 protocol.
