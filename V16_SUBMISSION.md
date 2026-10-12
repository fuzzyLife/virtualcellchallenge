# V16 published submission

V16 is published as `wG3GGbaTsl2dXFn7cQy1`, overall **0.141613857**, a gain of **+0.003730723** (2.71%) over V13 on the same panel and anchors. It is the best recorded project validation score. Rank was **562** in the user report and **561** at the API check; historical ranks are not directly comparable. See [V16_LEADERBOARD_RESULTS.md](V16_LEADERBOARD_RESULTS.md) for all components, raw comparisons and limitations.

- File: `experiments\v16_control_preserving\prediction.vcc`.
- Candidate: `residual10`; response: exact V13; generator: residual10.
- Size: 4,073,615,360 bytes (3.79 GiB).
- SHA-256: `aec15ac5e5ef38e835bb4206f62c431772a28690fe5823490f045c5f0319cc08`.
- Shape: 360,000 cells x 18,533 genes; 900 groups, exactly 400 cells each.
- Official prep and independent complete-payload/provenance/replay validation passed.
- Available public guards passed; complete six-metric public verification was unavailable.

Inspect the already published entry from Command Prompt:

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc status wG3GGbaTsl2dXFn7cQy1
```

The completed upload used model name `V13 response residual dispersion v16`. Repeating submit would create a duplicate experiment. Original generation-time upload instructions are archived in `experiments/v16_control_preserving/analysis/generation_submission_document.md`. [V16_RESULTS.md](V16_RESULTS.md) preserves the actual public tests and their limitations; the publication overlay does not alter their outcomes.
