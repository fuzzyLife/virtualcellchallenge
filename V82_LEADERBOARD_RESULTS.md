# V82 published: V66 recovery reproduced

The user-supplied CLI receipt reports **published** entry `sMObaQ7bhzac3RfNINSv`, model **V66 recovery after V69 - identical predictions**, overall **0.1508821875498449**, rank **550** at the time of the receipt.

All six reported scores and overall match V66 within 1e-15 on the same validation panel and anchors. The minuscule differences in overall and NMAE are floating-point precision, not evidence of a predictive change. This is the same prediction package, not a newly learned model. V66's old rank 518 and V82's rank 550 are different leaderboard snapshots; the score was reproduced while competitors' standings can change.

|Scaled score (higher is better)|V66|V82|Difference|
|---|---:|---:|---:|
|overall|0.15088218754984492|0.1508821875498449|-2E-17|
|pds|0.6517009514165119|0.6517009514165119|0E-16|
|mse|0.09657486167945828|0.09657486167945828|0E-17|
|nmae|0.08527061361861961|0.08527061361861947|-1.4E-16|
|fid|-0.02260859548432163|-0.02260859548432163|0E-17|
|reach|0.07962563498187135|0.07962563498187135|0E-17|
|jac|0.014729659086930079|0.014729659086930079|0E-18|

The recovery restores the best published project predictions after V69's 0.09526284481907588. It does not establish an improvement beyond V66, statistical significance, a new biological result or final-test performance. V70-V81 remain rejected under their original public criteria.

The local package SHA-256 is `be899434d81eac45bb8312077fe789088fe49b9b772971d589dcf5414e8f22f1`, 4,073,666,560 bytes. Prior preparation verified complete byte equality with V66 and official container validity. This recorder checks the saved evidence and current size; it does not rehash the large package or query the server. The receipt contains no server checksum or exact submission timestamp, so neither is invented.

The scoring job is `projects/vcc-prod/locations/us-west1/jobs/vcc-metrics-prod-smobaq7bhzac3rfninsv`. Upload duration was 965 seconds. Partition `val`, panel `vcc2026-val-1`, anchors `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`.

Reproduce with `python record_submission_result_82.py`, then `python record_submission_result_82.py --check`. Exact decimal strings and comparison CSV are retained in `analysis/v82_published`; publication input/output hashes are recorded there. Original preparation files keep their historical not-uploaded state; this publication record supersedes that status. No upload, API request, training or GitHub push is performed by the recorder.
