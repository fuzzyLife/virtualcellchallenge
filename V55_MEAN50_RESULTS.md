# V55: fixed half-blend diagnostic submission

This is the previously unsubmitted V37 mean50 response, produced in a new namespace. The user authorized a diagnostic after reviewing local-to-official validation mismatch. Original V37 failed its REACH guard and remains rejected in its original records. V24 remains the best published submission; V55 has no known leaderboard score.

## Completed paired tests

44 existing public populations,30 distinct context-target pairs,2 emission seeds,2 arms:176 predicted-population comparisons.40 original development comparisons were reused by checksums and136 new predictions generated. This is outcome-informed robustness testing on reused biology, not independent confirmation. D3 remains unopened. Rest D1 actual cells number200; other treatment populations400.

| Raw public metric | V24 | Fixed mean50 | Change |
|---|---:|---:|---:|
| pds | 0.775992063 | 0.775992063 | +0.000000000 |
| fidelity | 0.583803866 | 0.585123636 | +0.001319770 |
| reach | 0.218072360 | 0.216637676 | -0.001434685 |
| jaccard | 0.061525758 | 0.061428995 | -0.000096762 |
| nmae | 0.910677657 | 0.913649764 | +0.002972108 |
| common_expression_ratio | 37.450184307 | 37.335997886 | -0.114186421 |

Lower expression/NMAE errors are favorable; other raw metrics increase when favorable. These values are not anchored leaderboard scores. Equal targets within context, contexts within study and studies overall; repeated target populations are not independent competitors. Expression uses the existing panel jackknife calculation and retains its tiny-denominator limitation.

## Seed and case sensitivity

| Seed offset | FID change | Expression ratio relative change | REACH change |
|---|---:|---:|---:|
| 0 | +0.000964310 | -0.3063% | -0.006673848 |
| 550000 | +0.001675230 | -0.3035% | +0.003804479 |

Fidelity improves under every context-target removal: **True**.
Expression improves under every context-target removal: **True**.

The separate diagnostic technical bounds pass: **True**. These bounds permit an informative upload; they do not establish improvement or retrospectively relax the original V37 rule. Original strict REACH failure is retained. The full seed and30-case influence values are in assessment.json and independent_evaluation_audit.json.

## Study-specific trade-offs

| Study | FID change | REACH change | Expression ratio change |
|---|---:|---:|---:|
| CD4 | +0.005041417 | -0.000268464 | +0.000004998 |
| H1 | +0.000082805 | +0.010614948 | -0.000486223 |
| Orion | -0.001164911 | -0.014650538 | -0.342078039 |

CD4 drives most of the fidelity gain; Orion fidelity declines. Expression improvement is concentrated in Orion and its unstable HCT116 ratio denominator. These results do not establish gains across biological studies; the diagnostic status is retained. study_summary.json is a post-outcome descriptive summary and does not alter predictions or gates.

## What is packaged

For each context/target/output gene, effect = float32(V24 +0.5*(submitted V36 retrieval_contrast25 -V24)). This averages log2 effects, not integer matrices. No new network is trained. Unsupported/equal parent coordinates retain V24 exactly. Target-own effects remain excluded from these vectors and the unchanged generator applies its0.1 repression factor once.

The same V24 residual10 emitter, guide-stratified sampling, original control libraries, production seed20261011 and dual-moment calibration are used. The original V24 and V36 files remain unchanged. Synthetic cell IDs inherited from the validated generator use a v36 prefix; context/target labels and official axes are explicit.

Upload file: experiments/v55_mean50_diagnostic/prediction.vcc; **4,073,697,280 bytes**.
SHA-256: b877efdd31565ed7a58bcf8f26384a8837e546221fab5733ab0602e5db713c1b

Independent audits verify all528 saved-array DE metric values,20 profile-expression panels,30 target-removal expression calculations,5 unchanged null comparisons, the full production response blend and exact fallback. Preflight tests24 full-axis groups. Full output validates360,000 cells by18,533 genes,900 context-target groups with400 cells each. Every row library/provenance is checked. Official prep passes; full archive decompression compares every CSR value/index/pointer and gene/context/target label. Original control-null inflation persists; this experiment does not fix it.

## Submit from Windows Command Prompt

~~~bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v55_mean50_diagnostic\prediction.vcc -m "Fixed V24 V36 half-blend diagnostic v55" --wait
~~~

No automatic upload or Git push was performed. This is a diagnostic upload with uncertain score, not a recommended replacement for the V24 reference based on current evidence.

Reproduction: [V55_MEAN50_STRATEGY.md](V55_MEAN50_STRATEGY.md); root scientific driver run_v55_mean50.py, independent auditors audit_v55_evaluation.py/audit_v55_mean50.py, official packaging adapter package_v55_mean50.py. Source/array/receipt hashes are recorded in analysis/plan.json and the final delivery manifest. Large files remain local and regenerable.
