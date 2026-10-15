# V31 guide-resolved response teacher

This separately frozen pilot follows the V30 metadata audit: extracted CD4 donor2 resting/stimulated data contain separately assigned guide IDs. After requiring>=30 treated cells per guide and>=20 matched NTCs per lane,34 targets are represented per condition;22 resting and21 stimulated targets have multiple eligible guides. H1/Orion labels are predominantly co-delivered guide pairs and cannot identify each guide's separate coefficient.

## Fitting and calibration

Use only raw public CD4_D2Rest and CD4_D2Stim48 counts. Normalize each cell by its projected competition-axis count library. For each target/guide, calculate empirical expression means, sampling variances and lane-matched NTC means/variances. Source own-target suppression supplies an uncertain RNA knockdown proxy; the fixed prior is0.9, with correction gain bounded0.8-1.8. It is not an observed query dose or calibrated physical efficacy.

Pool guide responses with a random-effects variance that retains between-guide disagreement; cap precision ratios10:1. Confidence uses the normal approximation and total-cell reliability. Own genes are always absent from downstream response confidence. Compare with genuinely pooled all-guide raw means, guide labels shuffled within target/lane (preserving group sizes), and shuffled target fingerprints.

Calibrate the correction strength0.125/0.25/0.5 by withholding each entire actual guide identity across all lanes in D2 only. Baseline response estimates exclude all CD4 summaries during this calibration, avoiding donor-alias leakage. Shared source controls and the same culture make this within-study calibration, not independent validation. The held guide's downstream expression or suppression is never an input to its prediction. The fixed0.9 prior avoids learning a dose prior from held guide labels.

## Development boundary

Use previously spent H1 cases with H1 excluded and Orion cases with BOTH HCT116/HEK293T excluded. In CD4 donor1 resting tests, the baseline excludes every CD4 summary, while the new teacher uses raw donor2 counts. This is **donor transfer**, not whole-CD4 exclusion or unseen-lineage validation. Donor2 stimulated cells are training/calibration only, never validation. The remaining donor3 reserve stays untouched unless a separate count finalist passes its gates.

A query-control compatibility gate activates corrections only where CD4 controls are closer than non-CD4 source controls. Within-CD4 state weights use pairwise-valid rank correlations. These are source-selection hypotheses and do not identify anonymous A/B/C literally.

The first screen compares approximate mean CPM and clipped+/-3 smoothed0.1 log2FC against genuine held populations and independent reference controls, excluding own genes. It is a mean screen, **not**400-cell Wilcoxon testing or private anchor scaling. The same0.1 own-transcript factor and total normalization are applied to baseline/candidate mean profiles.

Frozen requirements:>=1% overall squared-LFC error improvement,>=0.5% absolute-LFC error improvement,>=5% D1 squared-error improvement, no family>2% worse, mean-CPM error<=1.01times baseline, and both error measures better than pooled, guide-shuffled and target-shuffled controls. A failed mean candidate does not advance to generated counts or reserve access.

## Reproduction and provenance

Final plan SHA256:cc2ea6ff449a8db040028f64e38b46fb2b2c7f7be5bfb7ac9226a095ef17b37c. Inputs, original extraction records, scientific code and runtime are hashed in the plan. V24/V30 and all previous frozen files are preserved.

```bat
"C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe" test_v31_guides.py
"C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe" run_v31_experiment.py prepare
"C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe" run_v31_experiment.py summary
"C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe" audit_v31_experiment.py
```

The already completed `freeze` stage is run only when creating a new untouched experiment. Source-only guide calibration precedes all development scores. No automatic upload is performed. Large NPZ arrays are regenerable and excluded from Git; small plans, labels, source row identities and verification results are versioned. Source identifiers are public cell/guide metadata, not private participant information.
