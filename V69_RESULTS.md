# V69 count-bridge release

Upload artifact: `experiments/v69_count_bridge/prediction.vcc`. Candidate: **hybrid25**. Official preparation and complete-payload validation passed. Public confirmation support: **False**. Overall release support including null checks: **False**. This file is **experimental**; no improved hosted score is claimed. V66 remains the highest published project score, 0.150882188.

## Change and decision

The V66 response means are preserved as model inputs. An independently implemented Gamma-Poisson count bridge changes raw integer counts through posterior intensities, binomial removal, Poisson addition and beta/gamma dispersion coupling. It is a working probabilistic model, not a claim that all raw cells are NB distributed. Source target/gene dispersion changes use split-view sign support, minimum molecule/cell counts and study-balanced shrinkage. All panel genes abstain from dispersion fitting. Shuffled evidence is an inference ablation, not a separately trained shuffled model. See [frozen protocol](V69_STRATEGY.md).

For hybrid25, exactly100 of400 cell positions receive bridge counts;300 retain the exact parent V66 rows. The mixture does not inherit the pure bridge's null identity. The mean_bridge/shape_bridge candidates use the full bridge population. No exact original library-size repair follows the bridge; final count validity and normalization are measured independently. Production samples only the selected100 bridge rows, using the same row-independent probability law with a documented new random-stream realization. Earlier code/metadata and an unused partial output are preserved; fitted arrays and selection rules are unchanged. Fresh preflight and deterministic replay verify the deployed sampler.

Development-eligible arms: []. Predeclared fallback used: **True**. Selection is locked before the remaining target cases are scored. Failed confirmation checks: profile, worst_study, fidelity, reach, jaccard, nmae, three_gains. Matched-shuffle profile comparison: **True** against **shuffled_hybrid25**. Unresolved metrics: ['nmae']. Confirmation cannot change the selected candidate. A separately frozen pre-development correction adds a shuffled quarter-mixture control, because the original full shuffled bridge was not a matched control for hybrid25. It changes no candidate or parameter. The confirmation assessor also preserves undefined DE measurements as unavailable and failed guards rather than attempting numeric comparisons with null.

Candidate paired equal-study profile error / V66: **1.107517** (+10.75%). These are reused-public raw results, not the leaderboard's anchored accuracies. Raw error and noise-corrected expression are both shown, because realistic sampling changes prediction variance.

## Confirmation

| Diagnostic | Fold-specific V66 | hybrid25 |
| --- | --- | --- |
| pds | 0.755556 | 0.822222 |
| fidelity | 0.584020 | 0.578927 |
| reach | 0.151328 | 0.137366 |
| jaccard | 0.068002 | 0.061966 |
| nmae | unavailable | unavailable |
| profile_square_error | 64.042135 | 71.095968 |

Higher PDS/FID/REACH/JAC is favorable; NMAE and raw squared profile error are lower-better errors. Missing metrics remain unavailable. Study profile ratios: {'H1': 1.0379688273352645, 'Orion': 1.171207171900939, 'CD4': 1.1133748661215939}. Different cases/target panel sizes prevent direct comparison to historical raw PDS numbers or hosted scores. Target confirmation reuses the existing public cohorts and is not independent-study validation.

## Mean versus decoder diagnostic

| Arm | Raw profile error | FID | REACH | JAC |
| --- | --- | --- | --- | --- |
| bridge_oracle | 138.313727 | 0.762463 | 0.991107 | 0.108748 |
| bridge_pred_oracle_shape | 140.797435 | 0.264708 | 0.399447 | 0.009475 |
| old_oracle | 92.772616 | 0.704057 | 0.944115 | 0.056258 |
| old_pred | 95.813455 | 0.520174 | 0.309471 | 0.040292 |
| real_fit | 92.581418 | 0.862441 | 0.970246 | 0.391633 |

Diagnostic representation checks: {'profile': False, 'fidelity': True}. Ten target/context pairs use distinct treatment fitting/scoring halves (100 or200 cells per synthetic prediction). Oracle means and shapes use the fitting half only; the scoring half is separate. The measured fitting-half mean is noisy and not a perfect population oracle. The bridge diagnostic additionally clips bulk ratios to the production +/-3 log2 range and cannot activate a gene with zero fitted control probability; the old-oracle input is not constrained identically. These limitations can contribute to reconstruction error and prevent attributing every difference to distribution family. Old-oracle generation receives both measured CPM/bulk means; bridge-oracle generation uses fitted bulk means and marginal shape. Realized means differ, so this is not a strictly orthogonal shape-only comparison. The real_fit arm is the separate treatment fitting half itself, not duplicated scoring cells. These labelled-data arms cannot be used for A/B/C inference or production selection. The source-trained deployable candidates use no query treatment labels. This diagnostic fails the raw-profile representation guard; that failure is retained. It does not justify claiming a generally superior decoder.

## Count-null checks and sampling correction

| Context | Representative template | V66 | Full bridge | Hybrid25 | Near-zero bridge | Original prefix bridge |
| --- | --- | --- | --- | --- | --- | --- |
| H1 | 0 | 1553 | 0 | 565 | 0 | 0 |
| HCT116 | 0 | 1963 | 0 | 678 | 0 | 578 |
| HEK293T | 0 | 1539 | 0 | 284 | 0 | 1041 |
| CD4_Rest_D1 | 0 | 1229 | 0 | 527 | 0 | 0 |
| CD4_Stim48_D2 | 0 | 2251 | 0 | 984 | 0 | 0 |

The original frozen null used the first400 query-fit controls. It can overrepresent source batches; exact count identity consequently preserved578/1041 DE calls in HCT116/HEK293T. That observation prompted a separately frozen [representative-null addendum](experiments/v69_count_bridge/analysis/representative_null_plan.json): fixed depth-stratified sampling of400 distinct query-fit controls, with no seed search. Original outcomes are retained. This addendum changes no fitting or candidate-selection rule. All arms in the additional comparison use the same selected rows. Zero calls in these comparisons do not establish a universal zero false-positive rate. The hybrid retains old-generator artifacts even when the full bridge passes.

## Noise-corrected expression

| Stage | Arm | H1 signed numerator/case | Orion | CD4 | Equal study |
| --- | --- | --- | --- | --- | --- |
| development | v66 | 0.000578 | 0.000508 | 0.006379 | 0.002489 |
| development | mean_bridge | 0.000719 | 0.001393 | 0.007514 | 0.003209 |
| development | shape_bridge | 0.000718 | 0.001393 | 0.007514 | 0.003208 |
| development | hybrid25 | 0.000628 | 0.000788 | 0.006723 | 0.002713 |
| confirmation | v66 | 0.001464 | 0.000395 | 0.003783 | 0.001881 |
| confirmation | hybrid25 | 0.001497 | 0.000549 | 0.003865 | 0.001970 |
| confirmation | shuffled_shape | 0.001578 | 0.001319 | 0.004046 | 0.002315 |
| confirmation | shuffled_hybrid25 | 0.001498 | 0.000554 | 0.003862 | 0.001972 |

These are signed, jackknife-corrected expression numerators averaged per case within context, then equally over contexts/studies. Lower is favorable; a negative estimate may arise from noise correction. Denominators/rho remain in the original panel receipts, and unresolved denominators are not floored into a claimed overall score. This descriptive analysis does not override the frozen raw-profile selection rule.

## Source evidence and limitations

| Raw source | Eligible target conditions | Supported target/gene coordinates |
| --- | --- | --- |
| H1 | 25 | 133950 |
| HCT116 | 234 | 690577 |
| HEK293T | 252 | 1069156 |
| CD4_Rest_D2 | 34 | 127839 |
| CD4_Stim48_D2 | 34 | 131730 |

HCT116/HEK293T share the Orion study; both CD4 conditions share the CD4 study. Entire recipient studies are excluded from both mean and dispersion evidence in public tests. A source view split uses alternating cells, not independent guides or biological replicates. The dispersion moments are not a Bayesian posterior or NB maximum-likelihood fit. Depth conditioning, clonal/batch effects and sparse counts can bias them. The new coupling preserves empirical control dependence only approximately after shape changes. No new biological dataset or untouched D3 reserve was opened. More seeds or reused public cases do not solve the validation-transfer gap observed with V66.

## Artifact and verification

- Rows/genes/groups: **360,000 /18,533 /900**,400 cells per target/context.
- Payload size: **3,992,545,280 bytes**. Nonzero entries: **2,297,638,080**.
- SHA-256: `faf909db096d55d21e2fc1d0f194ea18dc94e3a8c42acd1c57673be9366e06b4`.
- Independent payload audit: **360,000 rows**, all stored counts/indices/pointers, all control identities, **12** deterministic regenerated groups; **270,000** unchanged parent rows checked.
- Public result audit: **190** count populations and **1520** manually reproduced tied-rank coordinates, independent BH/LFC/profile checks.
- Maximum final library: **52,469**. Cells with changed total count: **89,034**. Library changes are allowed; exact original totals are not claimed.
- Six probability/count numerical tests and the quarter-population/cache equivalence checks pass. Package validity is not evidence of hidden accuracy.

The separate official environment runs unchanged `vcc prep`, following lossless float32/int64 staging. Every value is compared during staging, and the decompressed official payload is compared to the generated H5AD. Source arrays, model arrays and count populations are ignored regenerable outputs; compact receipts and root scripts are committed locally. No automatic upload or Git push occurs.

## Submit from Windows Command Prompt

```cmd
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v69_count_bridge\prediction.vcc -m "V66 probability-coupled count mixture v69 experimental" --wait
```

Keep the resulting entry ID. If waiting is interrupted after upload, use `vcc status ENTRY_ID --wait` rather than creating a duplicate submission. This release has no hosted result yet.
