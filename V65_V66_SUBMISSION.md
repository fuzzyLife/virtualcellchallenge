# Selected V66 submission: experimental test

Reviewed all 18 published submissions and tested four response-estimator variants across V65/V66. The fixed comparison selects **V66 `assay75`** under the frozen rule. **V63 (0.150422795) remains our best published reference.** This new file has no private leaderboard result. Public all-guard support: **False**. Failed checks: **reach**.

## Why this file

The trajectory favors pooled intervention evidence, calibrated coverage, context-matched donor responses and the recent assay-representation correction; ordinary blending/stacking did not improve on V24. V65/V66 isolate source normalization from correction strength while retaining the same count generator. The selection rule was committed before V66 count outcomes: prefer all-guard passers, otherwise only bounded diagnostic candidates with at least two raw improvements and no study profile degradation over 3%, then the smallest paired generated-profile error ratio. See release_selection.json for every finalist's outcome.

Selected count-profile error / V63: **0.996334814**. Raw metric changes: PDS +0.002976190; fidelity +0.001098608; reach -0.014431718; overlap -0.000033502; NMAE error -0.002445362. The generated-profile gain is concentrated in the held Orion study (HCT116/HEK293T); H1 and CD4 profile errors are slightly worse. These small reused-public-data trade-offs cannot forecast the official scaled overall. The 655-pair broader mean check is mixed and does not establish a clear improvement. Uploading measures this uncertainty; it is not a promised rank gain.

## Upload from Windows Command Prompt

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v66_assay_strength\prediction.vcc -m "Assay-strength response v66 experimental" --wait
```

Use this CMD syntax in the same terminal as previous submissions. If the wait is interrupted after upload, check the returned entry ID with `.\.competition_env\Scripts\python.exe -m vcc status ENTRY_ID --wait`; do not upload a duplicate. This release has not been uploaded automatically.

## Artifact verification

- File: `experiments/v66_assay_strength/prediction.vcc`, **4,073,666,560 bytes**.
- SHA-256: `be899434d81eac45bb8312077fe789088fe49b9b772971d589dcf5414e8f22f1`.
- Official local prep passes: 360,000 cells, 18,533 genes, 900 context/target groups, 400 cells each.
- Every packaged sparse count/index/pointer and all axes match the generated source.
- All original control identities, 400 distinct templates/group, guides and exact libraries checked; all counts are nonnegative integers within cell limits.
- All 900 production response vectors replay with every prior-margin constraint checked; 12 full production groups regenerate exactly.
- Ten additional primary public populations freshly replay profile/DE arrays and independent reach tie ordering; complete serialized audits cover 480 comparison records across both finalists.

These checks establish file/provenance integrity, not predictive accuracy. Scientific source and compact receipts are locally committed; large data/arrays/packages are intentionally ignored by Git. No Git push.

[All 18 published submissions](V65_SUBMISSION_HISTORY.md) · [V65 results](V65_RESULTS.md) · [V66 results](V66_RESULTS.md) · [V65 protocol](V65_STRATEGY.md) · [V66 protocol](V66_STRATEGY.md).
