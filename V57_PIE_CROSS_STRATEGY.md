# V57: cross-dataset PIE with genuine response evidence

This experiment implements the three steps authorized on 18 October 2026: genuine saved training evidence, all three retained heads, and source-fitted count calibration against V24 and shuffled controls. It does not promise a new leaderboard gain. Previous plans and packages remain unchanged.

The frozen machine-readable protocol is `experiments/v57_pie_cross/analysis/plan.json`. Eight root scripts are SHA256-sealed there. Their syntax was checked before the freeze. Any implementation repair must preserve the previous plan and explicitly record that no outcomes selected the repair.

## Evidence and conditioning

Use the published cross-dataset checkpoint `f6da9607226139a1040b6992736118c9b5c7c313`, SHA256 `e60a1304d074706405939bde94c5edf9851e78070f2d584c6b5025a002ed9a7c`, and the original saved split SHA256 `8d050f28aebf482135d42f758f6b440e0333b3dec8efc0dddab9dcc3099d4be4`. The source/code revisions match the V56 audit. No updated default asset reference substitutes for a saved checkpoint reference.

The collector retrieves 37 checksum-pinned files, about 17.12 GiB, into a resumable temporary cache. Existing V18 dense embeddings and full NCBI target token sequences are reused read-only. It does not download the entire 43 GB NCBI feature array or any new raw source expression matrix.

Stream all 64,260 original training rows across the four contributing datasets. The LFC clipping sample uses the exact author RNG, order, 2,048-row chunks and percentile. Global normalization uses every contributing gene and training row, accumulated in float64. The author evidence key must match the checkpoint. Only the 300 query targets' sufficient statistics are retained; this is explicitly a projection of the full saved evidence, not an assertion that a complete author cache was reconstructed. Independent arithmetic checks cover means, variances, flags, provenance and exclusions.

For H1 diagnostics, subtract that target's VCC25 source records. For HCT116/HEK293T diagnostics, subtract its matching line's records. This excludes direct query-response evidence but does not untrain weights already exposed to these source studies. CD4 has no matching context in these training metadata. Query contexts are presented as anonymous, with their real NTC mean log1p CP10k profile and missing context text. Context-response blocks are zero. The 154 seen targets and 146 unseen targets are reported separately.

All three outputs are retained: DE probability, conditional log2FC and mean per-cell log1p CP10k shift. Strict state-dictionary loading and repeated FP32 forwards check numerical integrity. Gene queries are restricted to measured, description-covered genes on contributing training axes; other genes retain the baseline rather than becoming zero-valued observations.

## Source fitting and transfer evaluation

Reuse the unchanged V52 disjoint cell-view row definitions. There are 366 eligible context-target cases: H1 12, HCT116 139, HEK293T 147, resting CD4 D2 34 and stimulated CD4 D2 34. These public populations are development data used in previous research, not newly blind biology. D3 remains unopened.

Fit calibration on H1/HCT116/HEK293T only. A whole-target hash split selects among three fixed ridge strengths; the selected model refits these source cases. No CD4 treatment values fit or select the decoder. Real PIE-head features, an independently trained target-shuffled control and an evidence-only ablation receive the same fitting procedure. The target shuffle preserves the seen/unseen strata. Train-constant features are disabled during normalization.

Three regression outputs use the model's three heads and donor statistics as features: a bounded residual to the V24 hold-study mean predictor, a rank effect, and a mean-log-expression shift. Rank effects are the original tie-aware probability shift, not a pseudocount transformation. The raw PIE mean-shift head is also assessed against a zero-shift baseline. Report source diagnostics and CD4 transfer separately, with seen/unseen results for both.

The V24 comparator is the existing public hold-study implementation, with the recipient study excluded from its response-source fit. This is not a re-upload or a claim to reproduce a private context-specific metric.

## Integer-count decoder tests

For each source context, choose two seen and two unseen metadata-eligible targets by a fixed hash. For each CD4 context choose four seen and four unseen targets. No outcome determines this selection. Source diagnostics use one fixed seed; CD4 comparisons use two seeds. This yields 176 four-arm population comparisons when complete.

The four arms are V24, calibrated PIE means with the unchanged V24 residual emitter, calibrated means with predicted rank-based dispersion, and the independently trained shuffled-head model with the same joint decoder. Real treatment populations contain up to 400 distinct available cells; generated populations contain 400 cells. Sample NTC templates and independent references without overlap. Seeds are repeated sampling, not independent biological replicates.

The joint decoder mixes compact Poisson counts and empirical control variation with a gene-specific nonnegative weight. Predicted mean changes and rank effects determine that weight through a Gaussian approximation. Unsupported or sign-inconsistent effects retain empirical variation. A zero intervention with zero rank effect anchors to empirical controls. Global mean/bulk moments and row library totals are enforced after the mixture. The approximation is tested end to end; it does not assume that mean and U determine tie-corrected significance.

Score actual integer counts after rounding and normalization: direction fidelity/reach, significance overlap, LFC error, pooled expression-profile MSE and delete-one corrected expression error. Evaluate raw DE confidence against real reference DE calls, excluding the perturbed coordinate. Run real-control null comparisons. Missing NMAE stays missing. The experiment has no private replicate anchors, so raw local scores are not converted into an invented official overall score.

## Runtime and commands

Python 3.12, the existing V20 scientific runtime plus original V18 PIE dependencies. Enable UTF-8 mode on Windows: otherwise the upstream default-text split reader misreads Greek drug identifiers. The supplied metadata is not renamed or altered to bypass that mismatch.

```bat
set PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe
%PY% -X utf8 collect_v57_pie_assets.py
%PY% -X utf8 build_v57_pie_evidence.py
%PY% -X utf8 prepare_v57_pie.py prepare
%PY% -X utf8 audit_v57_pie_cached.py evidence
%PY% -X utf8 infer_v57_pie.py --source H1
%PY% -X utf8 infer_v57_pie.py --source HCT116
%PY% -X utf8 infer_v57_pie.py --source HEK293T
%PY% -X utf8 infer_v57_pie.py --source CD4_Rest_D2
%PY% -X utf8 infer_v57_pie.py --source CD4_Stim48_D2
%PY% -X utf8 fit_v57_pie_decoder.py train
%PY% -X utf8 fit_v57_pie_decoder.py evaluate
%PY% -X utf8 evaluate_v57_pie_counts.py run
%PY% -X utf8 evaluate_v57_pie_counts.py summarize
%PY% -X utf8 audit_v57_pie_cached.py artifacts
```

The initial freeze command is `prepare_v57_pie.py freeze`; it refuses to replace an existing plan. Downloads and per-context inference resume with receipts. Large arrays stay in the temporary cache and can be regenerated; root scripts, compact result tables, provenance and checksums belong in Git. Do not commit weights or generated count matrices.

The independent auditor is invoked through `audit_v57_pie_cached.py`, which memoizes immutable JSON reads. The original frozen auditor, its formulas, samples and hash remain unchanged. This avoids repeatedly reading the same gene-description metadata for each gene. Reused source embeddings and author source are separately verified by `audit_v57_reused_assets.py`.
