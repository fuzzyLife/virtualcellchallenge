# V79 results: original development criteria fail; no new-model release

The protocol and executed code are recorded in`experiments/v79_count_innovations/analysis/plan.json` and`execution_protocol.json`. Ten reused public development target/context populations (first two distinct targets in each of five contexts) compare V66 with each genuine arm and its expression-bin gene/target-shuffled shape evidence control. Study means are balanced equally across H1, Orion and CD4. All raw profiles and discrimination are exactly V66; those invariants do not protect expression sampling corrections or DE means.

|Arm|Raw FID|Raw JAC|Raw REACH|Available-study raw NMAE|Fallback fraction|
|---|---:|---:|---:|---:|---:|
|v66|0.579503341|0.062351815|0.214770418|0.932783543|0.000|
|nb_query|0.393546377|0.058853602|0.269042274|0.932699062|0.667|
|turnover50|0.197468321|0.016167664|0.303752975|0.932906061|0.000|
|turnover80|0.353110102|0.027838835|0.250049192|0.932906084|0.000|
|source50|0.205248313|0.017201487|0.284904733|0.932860533|0.000|
|nb_query_shuffled|0.560173176|0.062444881|0.212314082|0.932800860|0.667|
|turnover50_shuffled|0.464059464|0.029411978|0.258782021|0.932688603|0.000|
|turnover80_shuffled|0.466400651|0.038000932|0.262160447|0.932626027|0.000|
|source50_shuffled|0.199519036|0.016306613|0.313618156|0.932924893|0.000|

No arm passes all original checks, so selection is null. No remaining-target confirmation, second seed, production population or upload file was generated for V79. These are unscaled public metrics, not the hosted anchor-scaled values. All public biology was reused; the tests do not establish a private leaderboard score.

The independent audit verifies all 90 saved baseline/candidate/control comparisons and 720 raw tied-rank coordinates, replays all four DE statistics and all BH arrays, recomputes PDS and verifies profile identity. Every fresh V66 baseline matched the frozen earlier cache. Fallback populations are original V66 and are not evidence of successful new shape calibration.

Query-fit controls were used to estimate shape; query-treated and held-out reference-control values did not fit the model.

Reproduce the saved audit using the pinned scientific environment:

```cmd
python -X utf8 audit_shape_trials.py --version 79 --stage development
```

The full per-arm Boolean checks and per-study values are preserved in`development_0_assessment.json`; the rejection is not a post hoc choice of a favorable metric. See`MARGIN_VARIANCE_STRATEGY.md` for mechanisms, original thresholds, count margins, dataset exclusions, runtime and limitations.
