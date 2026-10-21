# V90 diagnostic submission ready

**An upload-ready experimental file is available. Superiority over V66 is not established.** This is a post hoc leaderboard diagnostic using the unchanged, previously unsubmitted V86 interval25 configuration. V88 and V89 failed their original criteria and remain recorded as failures.

## File

`experiments/v90_interval_diagnostic/prediction.vcc`

- Bytes: 4,073,625,600.
- SHA256: `e042294a250ba151d533edb360c036d36522568240b4986438860ec000f78f81`.
- Cells: 360,000; genes: 18,533; contexts: 3; targets/context: 300; cells/target: 400.
- Official prep passed; the complete packaged sparse payload was compared independently.
- Every original selected template/library was verified; 12 full groups regenerate exactly.
- Uploaded: false; hidden score: unknown.

## Measured public tradeoffs

655 reused target/context pairs: desired-profile error ratio 0.996265839 relative to V66, raw PDS change -0.000536222. This is not an anchored Overall estimate. Strong shared/shuffled superiority criteria failed.

| Seed | Profile error reduction | PDS change | Fidelity change | Reach change | Jaccard change | NMAE change |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 0.2907% | +0.000000000 | +0.002504737 | +0.001932514 | +0.000892228 | -0.001844080 |
| 620000 | 0.3022% | +0.000000000 | +0.000325248 | +0.000000000 | -0.001626526 | -0.000782421 |

A descriptive target-level check found lower profile error on both donor seeds for 13 of 14 targets. This check followed candidate selection and did not alter production routing or predictions.

Both seeds improved profile error, fidelity and NMAE with unchanged discrimination. Jaccard worsened in the second seed, which failed V89's predefined overlap limit. The model does not have a validated overall improvement; this tradeoff is the reason for an explicitly diagnostic hosted test. No threshold was relaxed inside the original V88/V89 protocols.

The donor and broad public sets are extensively reused. Seed0 was known before this release choice; the technical seeds share biology. The original finite-iteration covariance fit and moment-based interval remain approximations. Production includes original V66 controls/templates/seeds and its unchanged residual10 count emitter.

## CMD upload command

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v90_interval_diagnostic\prediction.vcc -m "Conservative interval projection v90 diagnostic" --wait
```

Keep V66/V82 as the published reference until a new score is returned. A different rank on another date is not a controlled score comparison.

See V90_STRATEGY.md, the plan/production plan, full public/control comparisons in V88_RESULTS.md and V89_RESULTS.md, and the local official/independent validation receipts for reproduction. Large artifacts stay local and ignored; root code and compact evidence are versioned.
