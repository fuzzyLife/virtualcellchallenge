# V81 results: original development criteria fail; no new-model release

The protocol and executed code are recorded in`experiments/v81_residual_temperature/analysis/plan.json` and`execution_protocol.json`. Ten reused public development target/context populations (first two distinct targets in each of five contexts) compare V66 with each genuine arm and its opposite noise expansion control. Study means are balanced equally across H1, Orion and CD4. All raw profiles and discrimination are exactly V66; those invariants do not protect expression sampling corrections or DE means.

|Arm|Raw FID|Raw JAC|Raw REACH|Available-study raw NMAE|Fallback fraction|
|---|---:|---:|---:|---:|---:|
|v66|0.579503341|0.062351815|0.214770418|0.932783543|0.000|
|contract25|0.579694341|0.063036547|0.219481425|0.932794576|0.000|
|contract50|0.566792808|0.057294358|0.217165524|0.932713573|0.000|
|contract25_opposite|0.576115082|0.060917328|0.213724141|0.932692758|0.000|
|contract50_opposite|0.543138848|0.056674174|0.220947086|0.932685827|0.000|

No arm passes all original checks, so selection is null. No remaining-target confirmation, second seed, production population or upload file was generated for V81. These are unscaled public metrics, not the hosted anchor-scaled values. All public biology was reused; the tests do not establish a private leaderboard score.

The independent audit verifies all 50 saved baseline/candidate/control comparisons and 400 raw tied-rank coordinates, replays all four DE statistics and all BH arrays, recomputes PDS and verifies profile identity. Every fresh V66 baseline matched the frozen earlier cache. Fallback populations are original V66 and are not evidence of successful new shape calibration.

V81 additionally prespecified a control-only false-DE guard; it was conditional on development success and was not run. A lower-noise proposal is a calibration hypothesis, not a claim of biological realism.

Reproduce the saved audit using the pinned scientific environment:

```cmd
python -X utf8 audit_shape_trials.py --version 81 --stage development
```

The full per-arm Boolean checks and per-study values are preserved in`development_0_assessment.json`; the rejection is not a post hoc choice of a favorable metric. See`MARGIN_VARIANCE_STRATEGY.md` for mechanisms, original thresholds, count margins, dataset exclusions, runtime and limitations.
