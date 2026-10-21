# V80 results: original development criteria fail; no new-model release

The protocol and executed code are recorded in`experiments/v80_hurdle_controls/analysis/plan.json` and`execution_protocol.json`. Ten reused public development target/context populations (first two distinct targets in each of five contexts) compare V66 with each genuine arm and its expression-bin gene/target-shuffled shape evidence control. Study means are balanced equally across H1, Orion and CD4. All raw profiles and discrimination are exactly V66; those invariants do not protect expression sampling corrections or DE means.

|Arm|Raw FID|Raw JAC|Raw REACH|Available-study raw NMAE|Fallback fraction|
|---|---:|---:|---:|---:|---:|
|v66|0.579503341|0.062351815|0.214770418|0.932783543|0.000|
|hurdle_coupled|0.354683743|0.066451644|0.284659884|0.932751185|0.333|
|hurdle_independent|0.386247107|0.056977682|0.244339436|0.932756966|0.667|
|empirical|0.281209336|0.066913622|0.334173837|0.932656997|0.333|
|hurdle_coupled_shuffled|0.555433383|0.062674734|0.204728111|0.932762892|0.333|
|hurdle_independent_shuffled|0.553776036|0.062407164|0.205789682|0.932754127|0.667|
|empirical_shuffled|0.540601375|0.062961585|0.204573137|0.932759267|0.667|

No arm passes all original checks, so selection is null. No remaining-target confirmation, second seed, production population or upload file was generated for V80. These are unscaled public metrics, not the hosted anchor-scaled values. All public biology was reused; the tests do not establish a private leaderboard score.

The independent audit verifies all 70 saved baseline/candidate/control comparisons and 560 raw tied-rank coordinates, replays all four DE statistics and all BH arrays, recomputes PDS and verifies profile identity. Every fresh V66 baseline matched the frozen earlier cache. Fallback populations are original V66 and are not evidence of successful new shape calibration.

Query-fit controls were used to estimate shape; query-treated and held-out reference-control values did not fit the model.

Reproduce the saved audit using the pinned scientific environment:

```cmd
python -X utf8 audit_shape_trials.py --version 80 --stage development
```

The full per-arm Boolean checks and per-study values are preserved in`development_0_assessment.json`; the rejection is not a post hoc choice of a favorable metric. See`MARGIN_VARIANCE_STRATEGY.md` for mechanisms, original thresholds, count margins, dataset exclusions, runtime and limitations.
