# V24 context-matched donor response

V22 tested state-dependent response strengths and interactions; its development NMAE guard failed. V23 protected emitted means, passed development, then failed fidelity/Jaccard on new raw resting CD4 cells. V24 separately tests whether poor lineage/condition mean transfer is the limiting factor. Previous plans, numerical dependencies, model arrays and packages remain unchanged.

## Source and hypothesis

The existing global response pool uses weights K562: 2, CD4: 0.25, H1: 2, HCT116: 1 and HEK293T: 1, producing almost identical many-gene effects across contexts. Raw resting CD4 responses agree substantially better with author CD4 summaries than with the pooled predictor; the diagnostic finds no reversed-sign extraction bug. Private control A has higher pairwise-valid transcriptome rank correlation with CD4 than other sources; B/C do not. Compatibility is evidence for using a specialist, not proof of tissue identity or a guarantee of transferable causal effects.

Teacher inputs are the author's merged CD4 log2FC, SE and adjusted-P tables and two disjoint CE donor-pair tables already collected in reference/cd4_probe. Individual selected donors D1/D2 in the new raw release cannot be assumed to match or differ from the CE aliases; no raw-donor-disjoint claim is made. Source quality eligibility follows the existing build_cd4_model rule: >=100 cells, significant on-target knockdown, and no flagged off-target/neighbor/single-guide/low-expression issue. Source target knockdown calibrates effects to0.9, with gain clipped0.8..1.8 and LFC clipped+/-3. Target-own LFC remains zero for the established count emitter's separate0.1 knockdown factor.

Direction confidence is0.5 times the q<0.05 normal sign probability2*Phi(abs(LFC)/SE)-1. Agreement between both eligible donor pairs and the merged direction raises confidence to their minimum sign probability. Significant opposing donor directions abstain. Missing or invalid values abstain. No missing gene is a zero vote. Within-condition replicated contrasts are inverse-variance combined. Correction is evidence-weighted across conditions.

The additional null-equivalence arms use a bounded heuristic based on abs(LFC)+1.96*SE<0.1 before dose gain. This is not a calibrated equivalence posterior or formal TOST. The selected direction-only arm cdf50_de excludes these null corrections entirely.

## Rank routing and bounded correction

Use pairwise finite values, average tied ranks, genuine zeros and no pseudocount. Exclude all300 target-panel genes from routing. Compatibility gate is clip((bestCD4Rho-bestOtherSourceRho)/0.1,0,1), with other controls K562/H1/HCT116/HEK293T (excluding an evaluated held-out family). Within-CD4 condition weights are softmax(rho/0.1) on512 training-control genes with the largest rank changes. Held-out condition response rows are excluded, and the teacher's condition profiles during validation are Rest/Stim8hr only.

For each supported target/output gene:

```
correction = alpha * compatibility * confidence * clip(teacherLFC - baselineLFC, -0.5, 0.5)
newLFC = clip(baselineLFC + correction, -3, 3)
```

TeacherLFC receives the same basal gene/target expression gates as the baseline. Alpha is0.25 or0.5; selected cdf50_de uses0.5, giving maximum coordinate correction0.25. Own-target effects remain zero. Unsupported evidence or zero compatibility produces the exact baseline effect. Count generation retains the V16 residual10 emitter, dual-moment fitting, raw integers, original libraries and existing sampling seed20261011. This isolates the mean expert from count-emitter changes.

## Prospective experiment

The `cdf` candidate prefix labels the CD4 response expert. Frozen arms: v16, cdf25, cdf50, cdf50_de, their three alpha-matched target-shuffled controls, and an amplitude-only diagnostic. Target shuffles use fixed seeds24024/24025 and the same correction bounds; target-specific correction norms may differ, so this is not a perfect norm-matched control.

Development has28 previously studied genuine400-cell cases in H1/HCT116/HEK293T, plus eight200-cell raw resting D1 CD4 cases. The Rest cases overlap public summary biology and are **training-like calibration**, not blind donor validation. All CD4 Stim48hr response rows are removed from both baseline and expert fitting during development and confirmation. This baseline is the V16 mechanism refitted without those responses, not the published V16 package. All three non-CD4 families retain whole-family exclusion. HCT116 NMAE remains unavailable.

Release uses the predeclared available-metric guards, >=3 improvements including DE, worst-family fidelity/reach limits, and superiority to the candidate's alpha-matched target shuffle on >=2 of fidelity/reach/Jaccard. This mean-only experiment has no shape-improvement gate, specified before scoring. A single winner is selected by average metric-improvement ranks and locked before the fresh reserve. No re-selection after reserve is permitted.

Fresh reserve source: [author Zenodo21729034](https://doi.org/10.5281/zenodo.21729034), Marson2025_D2_Stim48hr_filtered.h5ad,3,474,151,389 bytes, MD5 24bcb324497d813a73fe85dfc559b8f5. Metadata-only power inspection found13 competition targets with>=400 cells. Freeze the first eight alphabetically with supported lanes,400 unique real cells each, and disjoint lane-stratified NTC halves. Full MD5, raw integer values, source IDs/labels, lanes, quality flags and128 complete original aligned rows are independently audited. Raw source data and big arrays remain outside Git.

The fresh raw condition's response rows are withheld from fitting, but it is the same public study and upstream differential-expression nuisance estimation may borrow other conditions. This is exploratory condition-transfer validation, not an independent-study test. The author's curated raw release contains999 targets plus controls; its selection bias is unknown.

## Production and reproducibility

A passing locked model is refitted on all three public CD4 conditions, including48hr, only after reserve. Its production baseline is the exact published V16 A/B/C effect arrays. Original challenge controls determine compatibility and condition weights; no hidden challenge treatment values enter. Identical effect groups reuse independently verified V16 count blocks losslessly. All900 groups retain400 unique, guide-stratified original templates; original row IDs/guides/libraries are saved. The package validator independently reconstructs every effect, checks all packaged sparse values/indices/pointers/axes and cell provenance, and replays12 full groups. Only the user uploads.

Exact scientific runtime: Python3.12.10, NumPy2.5.4, pandas3.0.6, SciPy1.18.1 and requirements-v23-science.txt. Pin BLAS to2 threads; pandas future.infer_string=False. Scientific scoring/generation/replay uses that environment. The unchanged competition environment performs official prep/validation/lossless decompression only. Existing environments were not modified.

Windows CMD (saved plans and receipts refuse overwrite):

```bat
cd /d L:\promec\Animesh\Download\Virtual
set "VCC_SCI_PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
"%VCC_SCI_PY%" test_v24_means.py
"%VCC_SCI_PY%" prepare_v24_experiment.py
"%VCC_SCI_PY%" run_v24_experiment.py development
"%VCC_SCI_PY%" assess_v24_experiment.py development
"%VCC_SCI_PY%" collect_v24_holdout.py
"%VCC_SCI_PY%" extract_v24_holdout.py
"%VCC_SCI_PY%" audit_v24_cd4_inputs.py
"%VCC_SCI_PY%" run_v24_experiment.py confirmation
"%VCC_SCI_PY%" assess_v24_experiment.py confirmation
"%VCC_SCI_PY%" audit_v24_experiment.py
```

The standard collector works from an empty data directory. The operational collect_v24_parallel.py recovery helper uses the saved four local partial prefixes; it is needed only to resume this interrupted download. Its transport amendment preserves source bytes and the original checksum; it changes neither targets nor model. For a passing final selection only:

```bat
"%VCC_SCI_PY%" fit_v24_submission.py
"%VCC_SCI_PY%" build_v24_submission.py --preflight
"%VCC_SCI_PY%" build_v24_submission.py
.\.competition_env\Scripts\python.exe package_v24_submission.py
"%VCC_SCI_PY%" validate_v24_package_recovery.py
"%VCC_SCI_PY%" document_v24_results.py
```

Plans, hashes, summaries, negative results and checks are versioned; large models/count caches/predictions are locally regenerated. A passing public experiment does not establish a higher hidden leaderboard score. The observed final decision and upload artifact are recorded in V24_RESULTS.md and V24_SUBMISSION.md.

The frozen fitter stored query gene labels as locally generated NumPy object strings. The original validator refuses object arrays with allow_pickle=False. The documented recovery wrapper verifies the creator/validator hashes and all three query file hashes, verifies every gene label, then permits object decoding only for those trusted local archives. It runs the original unchanged validator with every axis/count/formula/replay assertion intact. No numerical array or package is modified; the recovery receipt is included in the final audit.
