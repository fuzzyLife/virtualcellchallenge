# V15 release decision: retain V13

No new V15 package was created. All four deployable development candidates failed the frozen release criteria; see [V15_RESULTS.md](V15_RESULTS.md) for measured comparisons and [V15_STRATEGY.md](V15_STRATEGY.md) for reproduction. Reserved verification was not opened and no failed candidate was reselected.

The retained file is `experiments\v13_response_programs\prediction.vcc`, the existing V13 submission already published as `AvAsYkuxFSjbxYDCd0wJ` with overall 0.13788313368012328. Its SHA and official format were checked again after the V15 tests. **There is no reason to re-upload this identical package.**

- Size: 4,081,111,040 bytes (3.80 GiB).
- SHA-256: `2454a8b8ce694217962b89f1bac4ec0fb1bbd01a396c1ffe5cd780f4baa11591`.
- Shape: 360,000 cells x 18,533 genes; 900 groups of 400 cells.
- Artifact check: `experiments/v15_weighted_transport/analysis/baseline_package_check.json`.

If this file is needed for a different permitted evaluation, its Command Prompt command is:

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v13_response_programs\prediction.vcc -m "Assay-aware knockdown transport v13" --wait
```

No upload was performed by V15. The retained V13 model was exploratory when built; its published score is the observed evidence, not a guarantee on the final test set. A format check does not validate biological accuracy.
