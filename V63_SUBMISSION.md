# V63 experimental test submission

V63 is ready for a diagnostic upload. Public development gates failed, so this is **not a validated replacement for V24**. The source/assay correction improves expression-error diagnostics, with a remaining target-discrimination trade-off. Read [all tests and limitations](V63_RESULTS.md).

File: `experiments/v63_identity_guard/prediction.vcc` (4,073,656,320 bytes).
SHA256: `a075397d8c23fa2e63e1e99220a04bf6fcd500eaa1e94f8a0fc8eb445271048b`.
Official prep and the full package/count/axis/provenance audit passed, including 12 exact group replays. No hosted V63 score is known.

Windows Command Prompt:

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v63_identity_guard\prediction.vcc -m "Identity-guarded assay response v63 experimental" --wait
```

If waiting is interrupted after upload, use vcc status with the returned entry ID. Do not start a second upload. No upload or Git push was performed here.
