# V61 experimental submission

The fixed primary failed the recorded public support gates. This release is experimental and is not recommended as a validated replacement for V24.

Ready file: `experiments/v61_pooled_state_response/prediction.vcc`.

- Size: 4,068,689,920 bytes.
- SHA256: `cf4807ba59156619ae1bc2a315a87450bab5a1021c5c8cc569703b059118ad53`.
- Shape: 360,000 cells × 18,533 genes; 300 targets × 400 cells × A/B/C.
- Official VCC prep: passed.
- Complete packaged value/index/pointer and axis comparison: passed (2,365,303,864 nonzeros).
- Original control provenance/library checks: all 360,000 cells.
- Exact count replay: 12 full groups.
- Private leaderboard result: not submitted; unknown.

Windows CMD:

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v61_pooled_state_response\prediction.vcc -m "Pooled state-response transfer v61 experimental" --wait
```

This uploads an experimental diagnostic candidate. V24 remains available at `experiments/v24_cd4_response/prediction.vcc`; the new release does not overwrite it. [Measured tests](V61_RESULTS.md) and [reproduction](V61_POOLED_STATE_STRATEGY.md).
