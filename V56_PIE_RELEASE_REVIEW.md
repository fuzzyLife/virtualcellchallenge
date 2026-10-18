# V56: PIE release audit and the remaining gap to leading submissions

Audit date: **18 October 2026**. This is an investigation of the cloned release, saved model metadata, our V18 implementation and published project results. It does not train a model, download new model weights, produce a submission or establish a leaderboard gain. V24 remains our best published submission, overall **0.142997855**; V55 scored **0.138767074**.

## Findings that change the next experiment

**We have not fairly exhausted PIE.** V18 rejected a particular within-dataset checkpoint and two ways of blending its LFC predictions into V16 counts. It did not test the broader cross-dataset checkpoint, its actual response-evidence pathway, or its mean-shift head in a count decoder. The narrow experiment was useful and correctly rejected a harmful blend; it is not evidence that every PIE configuration fails.

The newly cloned repository is the same release already archived in V18. Local HEAD and the current GitHub API snapshot both identify `21a526add9ff3a81b4f962d9b16bbae570052c63`. All 162 corresponding archived files match after normalizing checkout CRLF to LF; 151 differ only in those line endings. Neither Hugging Face model repository revision has changed since the prior audit:

- Within dataset: `055c7a2cabdaff5659121c64fa2b609c231fc5be`.
- Cross dataset: `f6da9607226139a1040b6992736118c9b5c7c313`.

Thus the opportunity is an **untested configuration of existing weights**, rather than newly improved weights. Upstream fixes for saved data paths and source aliases are already in the archived revision. New default experiment asset references must not silently replace the older references saved in a checkpoint.

Sources: [Arc announcement, October 5](https://arcinstitute.org/news/pie), [author repository at the audited commit](https://github.com/ArcInstitute/pie/tree/21a526add9ff3a81b4f962d9b16bbae570052c63), [within-dataset model card](https://huggingface.co/arcinstitute/PIE_replogle_wdataset), [cross-dataset model card](https://huggingface.co/arcinstitute/PIE_replogle_xdataset). The pasted release agrees with the official announcement. The linked [preprint](https://doi.org/10.64898/2026.10.02.756297) was inaccessible through the browsing tool; no full-paper review is claimed here. Implementation conclusions below are checked against the author source.

## Actual target and gene coverage

The new audit downloads the exact saved training splits and asserts that their SHA256 values match each checkpoint's saved statistics. Coverage uses canonical gene symbols, rather than an assumption based on dataset names or an unverified alias conversion.

| Property | V18 checkpoint: within-dataset, K562 held out | Untested cross-dataset checkpoint |
|---|---:|---:|
| Training context-perturbation pairs | 2,603 | 64,260 |
| Our 300 targets present in the actual training split | **0** | **154** |
| Our output genes on contributing dataset axes | **6,290 / 18,533** | **18,260 / 18,533** |
| Response evidence for our targets | None | Potentially available for those 154 targets |
| Training sources | HepG2, Jurkat, RPE1 | Tahoe, Jiang, VCC25, Orion |
| Official 2026 competition score established by these model cards | None | None |

Cross-dataset output coverage is the union of its four contributing dataset gene axes intersected with our axis and the available gene descriptions. It is **not** a claim that every coordinate has significant supervised examples or predicts accurately. The union contains 31,675 genes. Its 18,260 competition coordinates leave 273 requiring a fallback or a separately validated extrapolation.

Actual training-target overlap by source is Orion **148**, VCC25 **13**, Jiang **6**, Tahoe **0**, with **154 unique** targets after removing overlap. Metadata vocabulary overlap is larger for VCC25/Jiang (16/7) but includes targets outside the training split. Tahoe perturbations here are drugs; zero exact gene-symbol overlap does not imply that chemical training is useless.

The cross-dataset checkpoint includes HCT116/HEK293T in Orion and H1 in VCC25. Those reused public studies cannot establish unseen-study performance for it. Jiang's saved contexts are stimulated A549/BxPC3/HAP1/HT29/K562/MCF7 lines. No CD4 context is listed in the four training metadata files; CD4 is therefore a useful potential transfer test, subject to accession/cell-provenance checks. An unopened donor is not automatically an independent study, and repeated public tests remain development evidence.

The within-dataset `k562/` directory denotes the **held-out** line. It is not a K562-trained checkpoint. The published four folds are independent models, not a reported ensemble. The model cards' AUPRC/direction/Jaccard results use their own benchmark definitions; they are not the six scaled VCC2026 scores.

## What our earlier adapter actually tested

The root `infer_v18_pie.py` uses the strict author state dictionary and original saved split. `serve_evidence` returns zero donor blocks for these unseen target identifiers and zero context-response blocks for a novel dataset/context. Our zero blocks were appropriate for that restricted configuration; they were not a substitute for evidence available to the cross-dataset checkpoint. Query control-expression features and available biological source tokens were still present.

We deliberately left anonymous A/B/C context text missing. Guessing a Cellosaurus identity from VIM or other markers is not a reliable replacement. This missing modality is allowed by the model but makes our use different from benchmarks with known cell-line descriptions.

The old code computes three outputs but persists only `lfc` and `p_de` in its inference NPZ files. `delta_p` is discarded. V18 then tests a 25% LFC blend and a confidence-gated blend on the covered genes, retaining V16 elsewhere.

V18's 28 public target-seed cases produced these raw metrics; these are not official scaled leaderboard scores:

| Metric | V16 baseline | PIE confidence blend | PIE 25% blend |
|---|---:|---:|---:|
| PDS | 0.77778 | 0.77778 | 0.65278 |
| Direction fidelity | 0.60648 | 0.60485 | 0.54239 |
| Direction reach | 0.50343 | 0.52895 | 0.38325 |
| DE Jaccard | 0.02237 | 0.02172 | 0.01862 |
| NMAE error, lower is better | 0.84953 | 0.86310 | 0.90355 |
| Relative expression worsening versus baseline | 0 | +6.61% | +650.30% |

The confidence arm failed the fixed acceptance rules and did not establish useful complementarity versus amplitude and shuffled-target controls. V18 was never submitted. Its plans, results and failure remain unchanged.

## The important distinction between PIE's heads

The author [loss implementation](https://github.com/ArcInstitute/pie/blob/21a526add9ff3a81b4f962d9b16bbae570052c63/src/pie/model/losses.py) trains separate objectives:

1. DE classification over tested genes, with class weighting.
2. Log2 fold-change regression **only on DE genes with valid, positive linear fold changes**; a separate direction penalty uses the same mask.
3. A binned mean-expression-shift objective, including the available non-DE coordinates.

A fold-change head trained conditionally on DE is not automatically a calibrated predictor of small or zero changes throughout the transcriptome. The unweighted all-gene interpretation is especially questionable for the unconditional blend. This is a plausible contributor to its failure, not an isolated causal explanation of the measured loss.

`p_de` is a tempered, class-weighted model output, **not a Wilcoxon p-value or BH q-value**. Applying a fixed probability threshold does not reproduce the competition's 400-cell significance calls. The old probability gate is an experimental blend rule, not a validated conversion.

`delta_p` means a change in **mean per-cell natural-log1p CP10k expression**. It is not detection probability, rank shift or Mann-Whitney U. There is no hidden standardization of the LFC output in `readout`: LFC is passed through unchanged. Our earlier code did not simply miss an output de-standardization step.

These scales also cannot be interchanged through one exponentiation:

```text
mean(log1p(normalized cell counts)) != log1p(mean(normalized cell counts))
```

Our scorer-style expression profile uses log1p of pooled counts normalized to 50,000. A useful decoder must check that profile as well as per-cell mean shifts, fold changes, ties and DE calls. PIE supplies population statistics, not the final 400 integer-count cells required for each competition group.

## Why our leaderboard score is far behind

The following comparison uses the user's latest pasted leaderboard, with the leader's components rounded to three decimals. It is not a new live ranking snapshot. All columns below are **scaled accuracy scores: higher is better**, including MSE and NMAE.

| Component | Our best, V24 | Illumina AI / iCell, user snapshot |
|---|---:|---:|
| PDS | 0.644 | 0.878 |
| Expression accuracy | **0.048** | **0.786** |
| LFC accuracy | **0.082** | **0.438** |
| Direction fidelity | -0.016 | 0.172 |
| Direction reach | **0.085** | **0.463** |
| DE overlap | 0.015 | 0.088 |
| Overall | **0.1430** | **0.4709** |

Expression, LFC accuracy and direction reach explain approximately **75% of the component-score gap**. A better perturbation ordering alone cannot close it. This is score arithmetic, not proof of which methods the leader uses. Searches did not locate a verified technical recipe tying iCell's 2026 leaderboard result to PIE. The PIE announcement likewise does not establish that the leading teams use it.

Our own diagnostics identify concrete weaknesses: V24's intended CD4 directions agree with held public DE directions only around 53-57%; synthetic zero-effect populations call 355-1,415 DE genes versus 0-1 in the corresponding real-null tests; and the reused local expression-ratio aggregate is dominated by one unstable denominator. V55's local fidelity improvement failed to transfer privately. These are documented in [V51](V51_PLATEAU_DIAGNOSIS.md), [V54](V54_METHOD_PERFORMANCE.md) and [V55's published result](V55_LEADERBOARD_RESULTS.md).

We have already trained mean/rank networks, tested population models and tried ensembles. The gap is not evidence that we forgot neural networks or majority voting. The untested opportunity is a materially broader pretrained response predictor with appropriate conditioning, followed by an independently checked conversion to counts.

## Recommended next experiment

This is a proposed experiment, not an executed performance result or a promise to beat V24.

1. **Cross-dataset PIE with its saved contracts.** Use the pinned checkpoint, saved source/alias references and hash-matched training split. Build genuine donor evidence from that split, with the query response excluded. Keep A/B/C context descriptions missing. Separate the 154 seen targets from the 146 unseen targets in every report.
2. **Retain all three outputs.** Test mean-shift predictions before synthesizing cells. Compare against V24, an evidence-only transfer baseline and a source/shuffled-target control. Evaluate both the saved LFC predictions on their intended DE domain and the mean head across measured genes. The first pilot should answer whether this checkpoint adds real response information.
3. **Separate overlapping-source diagnostics from transfer tests.** Orion/H1 tests can check numerical integration but cannot prove zero-shot performance. Use provenance-checked CD4 populations as a distinct setting. Any tuning on those public populations must be reported as development. Freeze choices before opening the reserved donor; include controls and multiple fixed generation seeds without counting seeds as new biology.
4. **Test a joint count decoder rather than another fixed LFC average.** Calibrate the relationship between predicted mean shifts, conditional LFC, DE confidence and integer-cell distributions using training sources only. Keep real query control distributions as the zero-response anchor. Check pooled expression profiles, sample-size sensitivity, ties, DE-set precision, direction and null calibration after row normalization and rounding. Uniformly restoring control variance already failed in V29; variation must not be changed indiscriminately.
5. **Gate any package on end-to-end evidence.** Require an improvement that survives target/context exclusions and comparison with simpler controls. If only the overlapping training studies improve, report that limitation instead of calling it a generalization result. Start with batch-one/chunked inference and existing cached features; the roughly 0.97 GB checkpoint does not require retraining the published eight-H100 configuration.

The immediate priority is **a correctly conditioned cross-dataset PIE pilot using its mean head**, before allocating another leaderboard submission to a small blend of the existing response vectors.

## Reproduce this audit

Root script: `audit_v56_pie_release.py`. Artifacts: `analysis/v56_pie_release/` contains ten bounded metadata snapshots with URL/size/SHA256 receipts, exact training-target coverage, code/model revision comparison and the rounded leaderboard-gap calculation. The largest new snapshot is the 1.64 MB training split. No author code is executed by this audit.

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe audit_v56_pie_release.py
```

Requires `requests`, Git, the unchanged local `pie` clone and the prior V18 metadata/source archive. To create the metadata snapshot in a fresh location where it does not yet exist, use `--fetch`; existing snapshots are never silently refreshed. Offline replay verifies the saved source hashes, recomputes target overlaps, asserts split hashes against checkpoint statistics and compares the archived author files with the clone. The user's cloned repository and saved HTML assets remain untouched. No old predictions, frozen plans or submission receipts are modified.
