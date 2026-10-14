# V24 submission status

Passed the locked raw-condition reserve.

Frozen candidate: cdf50_de. Upload-ready, independently verified; not uploaded.

File: `experiments/v24_cd4_response/prediction.vcc` (4,073,656,320 bytes).
SHA256: `42cd3906e17aa2e8e546efea44f1e726a38e577605f211ad96f1564362a20cc6`.

Windows Command Prompt:

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v24_cd4_response\prediction.vcc -m "Context-matched donor response v24" --wait
```

V16 remains the published champion until the host scores this file. Local tests do not guarantee a higher leaderboard result. If upload finishes and waiting is interrupted, use vcc status <entry-id> --wait instead of creating a duplicate submission.

Methods and exact reproduction: V24_STRATEGY.md. Actual tests/limitations: V24_RESULTS.md. Code is in the project root; small evidence is in experiments/v24_cd4_response/analysis. Large sources, arrays and packages are ignored in Git.
