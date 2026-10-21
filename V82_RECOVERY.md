# V82 recovery: exactly the V66 predictions

This release restores the strongest published project model after V69's regression. **It contains no newly learned predictions.** V70-V81 all failed their original public requirements; their results and strict rejections are retained. See `POST_V69_RESULTS.md` for the evidence and `MARGIN_VARIANCE_STRATEGY.md` for reproducible protocols.

File: `experiments/v82_v66_recovery/prediction.vcc`
Bytes: 4,073,666,560
SHA-256: `be899434d81eac45bb8312077fe789088fe49b9b772971d589dcf5414e8f22f1`
Status: byte-identical to V66; official container validation passed; not uploaded.

The original V66 entry `x8V3fUSgFIoAzHbNx4eX` scored **0.15088218754984492**. This is a previously observed score, not a new result. Original V66 public guards also had failures; those records are preserved. Equality of the complete package ties the recovery to the earlier full-payload audit of all 360,000 cells and 900 groups. Current rank and final-test performance are not guaranteed.

To restore those predictions as the latest leaderboard entry, use Windows Command Prompt:

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v82_v66_recovery\prediction.vcc -m "V66 recovery after V69 - identical predictions" --wait
```

This is a repeat of the V66 predictions and consumes a submission allowance. There is no scientific reason to re-upload it if you do not need to restore the current leaderboard entry. Do not interpret the V82 directory name as a new model experiment.

`prepare_v82_recovery.py` verifies that all twelve pilots completed and failed release, copies the unchanged champion, checks source/copy SHA-256, and runs the official container validator. It refuses to overwrite any existing differing output. Use the official CLI environment for this helper. `summarize_post69_trials.py` regenerates the complete comparison report. Code and compact evidence are committed locally; large arrays and upload packages are ignored by Git. No GitHub push or submission was performed.
