# V16 upload-ready exploratory prediction

Official prep and complete-payload validation passed. No leaderboard result is known. See [V16_RESULTS.md](V16_RESULTS.md) for actual public gains, missing metrics and limitations.

- File: `experiments\v16_control_preserving\prediction.vcc`.
- Candidate: `residual10`; response `v13`; generator `residual10`.
- Size: 4,073,615,360 bytes (3.79 GiB).
- SHA-256: `aec15ac5e5ef38e835bb4206f62c431772a28690fe5823490f045c5f0319cc08`.
- Shape: 360,000 cells x 18,533 genes; 900 groups, exactly 400 cells each.
- Available public guard verification: True; complete six-metric verification: False.

Run in Command Prompt:

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v16_control_preserving\prediction.vcc -m "V13 response residual dispersion v16" --wait
```

After an interruption following upload, use the printed entry ID with `vcc status ENTRY_ID --wait` rather than submitting again. No upload was performed here. V13 remains the best published reference until a new score is available.
