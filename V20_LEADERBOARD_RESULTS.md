# V20 published leaderboard result

V20 is published as `B92qai3cUvEHeobzp8Fk`, overall **0.139011567**, a change of **-0.002602290** (-1.84%) from V16 on the same panel and anchors. **V16 remains the best recorded submission.** Rank was **552** in the user report and **554** at the API check; historical ranks cannot establish improvement when the leaderboard changes.

## Scaled comparison

| metric | v16 | v20 | change_from_v16 | overall_contribution_change |
| --- | --- | --- | --- | --- |
| overall | 0.141614 | 0.139012 | -0.002602 |  |
| pds | 0.641553 | 0.639288 | -0.002265 | -0.000377 |
| mse | 0.046085 | 0.047823 | 0.001738 | 0.000290 |
| nmae | 0.080345 | 0.079577 | -0.000767 | -0.000128 |
| fid | -0.017296 | -0.016123 | 0.001173 | 0.000195 |
| reach | 0.084381 | 0.068935 | -0.015446 | -0.002574 |
| jac | 0.014615 | 0.014569 | -0.000046 | -0.000008 |

All six displayed components are scaled scores: **higher is better**, including MSE and NMAE. Expression accuracy and direction fidelity improve; PDS, fold-change accuracy, reach and significance overlap decline. Reach contributes **-0.002574297** to the overall change, **98.9%** of the net decrease. This is arithmetic decomposition, not causal attribution.

## Raw comparison

| metric | v16 | v20 | change_from_v16 | relative_change | higher_is_better | improved |
| --- | --- | --- | --- | --- | --- | --- |
| pds_cosine | 0.790561 | 0.789617 | -0.000944 | -0.001194 | True | False |
| expr_mse_unbiased_capped_norm | 0.961835 | 0.959587 | -0.002249 | -0.002338 | False | True |
| de_wilcoxon_lfc_nmae | 0.952086 | 0.952407 | 0.000321 | 0.000337 | False | False |
| de_wilcoxon_direction_fidelity_yield_raw | 0.507558 | 0.507932 | 0.000373 | 0.000736 | True | True |
| de_wilcoxon_direction_reach_raw | 0.153714 | 0.139739 | -0.013976 | -0.090920 | True | False |
| de_wilcoxon_sig_jaccard | 0.035660 | 0.035645 | -0.000015 | -0.000419 | True | False |

Raw DE reach falls **9.09%**. Raw expression error decreases slightly and raw fidelity improves slightly; raw NMAE worsens. The public reserve predicted reach improving by +.022051, but the hosted raw reach decreases by -0.013976. Passing the available public guards did not establish transfer to anonymous A/B/C. Missing HCT116 NMAE and unaudited pretraining overlap remain limitations; they must not be retroactively described as complete or independent biological validation.

## Interpretation and next test

The result supports retaining V16. It does not support another blind increase in ensemble weight. A separate pre-emission audit finds total absolute eligible LFC mass changes of -0.061% in A, +0.786% in B and +1.519% in C. Thus the mean correction did not uniformly shrink response amplitude. These sums include non-DE genes and provide no hidden truth: they cannot tell which genuine responses were lost. Gene-specific calibration, source/context mismatch and nonlinear count/significance effects remain hypotheses, rather than identified causes.

Before another submission, prioritize a fresh context/assay or donor holdout, verify pretrained-study overlap where possible, and diagnose response gains/losses per context and target on accessible treatment truth. Keep V16 generation fixed while testing a narrowly defined correction against matched magnitude/shuffled controls. Use a new reserved set for any new selection; the V20 reserve has been spent. No new candidate is built by this publication recorder.

## Provenance

The official API confirms publication, score, panel and anchors. It does not provide an uploaded-file checksum; the user's CLI path ties the audited local package to the entry. Original models, predictions, validation receipts and the generation-time `submission_ready.json` remain unchanged; that sealed receipt still records uploaded=false at generation. Current publication state is in `leaderboard_result.json`. Original report documents are archived byte-for-byte and mapped by `publication_manifest.json`. After rerunning the generation report writer, rerun this recorder to restore the publication overlay.

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc status B92qai3cUvEHeobzp8Fk --wait
.\.competition_env\Scripts\python.exe record_v20_result.py
```
