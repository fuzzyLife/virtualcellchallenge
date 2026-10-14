# V20 submission: independently validated and ready for manual upload

Artifact: `experiments/v20_stack_transfer/prediction.vcc` (3.79 GiB; 4,073,574,400 bytes).

SHA-256: `8370080f290febfbec59e25a54b175317f1a22609ee3c90faa1997774611ac50`.

The file has **not been uploaded** and has no leaderboard score. V16 remains the best published submission until the new result is reported. V20 is the bounded Stack mean ensemble that passed both the development and once-only reserved available-metric guards. The active V21 variation mechanism was rejected and is not included. Missing HCT116 NMAE, reused public biology and unaudited pretrained-study overlap limit the validation claim; a leaderboard improvement is not guaranteed.

## Upload from Windows CMD

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v20_stack_transfer\prediction.vcc -m "Stack-calibrated response ensemble v20" --wait
```

If waiting is interrupted after upload, use the entry ID returned by the CLI with `.\.competition_env\Scripts\python.exe -m vcc status <entry-id> --wait`, rather than starting another upload.

## What passed

- 644 full-count development comparisons across four research branches, followed by 80 reserved comparisons of the selected mean model and fixed controls. These are paired algorithm comparisons, not 724 independent biological replicates.
- All 900 frozen public source-choice records reproduced before private preparation. Every production query and all 587 source prompts were independently compared with original raw counts, labels, batches, sampling seeds, gene axes and missing padding.
- All 900 production correction vectors independently reconstructed from original V16 effects and generated expert arrays; own targets, unsupported responses and the .125 log2 bound verified.
- Official unchanged VCC prep and package validation. Complete packaged counts, indices, pointers and axes compared with the original generated H5AD: **2,368,754,215 nonzero values** checked.
- All **360,000 cells**, **900 groups**, 18,533 gene positions, original source-cell identities, balanced NTC guides and original integer library sizes verified. Twelve groups exactly regenerated, including the largest changed responses. 36 bitwise-unchanged V16 groups reused after verifying the original artifact.
- Safe plain-string metadata repair and pandas metadata compatibility recorded; all model parameters, numeric inputs and reserved decisions remain fixed.

Methods and reproduction: [V20_STRATEGY.md](V20_STRATEGY.md). Public results and limitations: [V20_RESULTS.md](V20_RESULTS.md). Release evidence: `analysis/submission_ready.json`, `independent_validation.json`, `production_input_independent_audit.json`, `reference_model.json`, `production_plan.json`, `production_response_changes.csv` and runtime/repair receipts under the experiment directory. Large raw data, weights, regenerable numeric caches and full sampling provenance are retained locally with hashes and excluded from Git.
