# V13: response transfer experiments

The goal is a separately versioned upload artifact after testing new response ideas against V11. V11 remains the best published entry (overall 0.137003974). Local source-panel scores cannot predict the hidden ABC leaderboard.

## What the research changes

Arc's October 1 data description confirms matched guides, time point and assay across contexts, strong target knockdown, and curated 400-cell groups. The response-decomposition preprint distinguishes generic stress/template responses from perturbation-specific and context-interaction components; baseline expression does not resolve all interactions. A separate tabular-model benchmark supports evaluating informative response features and simple models before increasing architecture size. Primary links and specific actions are recorded in `experiments/v13_response_programs/analysis/literature_review.json`.

We downloaded and verified the author's 95,350,546-byte RPE1 raw pseudobulk release. Its intervention symbols have zero overlap with the 300 competition targets after explicitly parsing the composite construct/promoter/transcript IDs. It shares 8,260 output genes, but those unrelated interventions are not direct evidence for our targets. The complete identifier audit is `audit_v13_rpe1.py`, including five published missing-Ensembl construct labels; the original coverage record is preserved. The reproducible downloader is `collect_v13_rpe1.py`; published MD5, SHA-256, and coverage are recorded in the analysis folder. RPE1 could still support future functional-prior research.

## Separate response experiments

1. **Response programs:** derive a 64-dimensional SVD basis from training-family V11 response estimates, whitened by median source uncertainty. Exclude reserved confirmation targets from basis fitting. Project onto 8, 32 or 64 programs and blend 25% or 50% with the original estimate. Keep missing gene coordinates masked and the target's own effect zero. This is an independently implemented test inspired by the decomposition paper, not a reproduction of its reported model.
2. **Robust source consensus:** use the median across measured lineages to limit isolated source effects when at least three independent families are present. Test LFC-radius floors 0.25 and 0.50, retaining original source disagreement in uncertainty. Separately test five Student-t reliability iterations. Correlated CD4 conditions remain one family.
3. **Assay-aware knockdown quality:** an audit found 149/291 HEK293T contrasts below 50% estimated knockdown, while all 41 zero-KD K562 rows have an unmeasured own-gene readout. The old dose correction amplified these missing measurements as though they were weak measured knockdowns. Replace unknown/unquantifiable row KD with a neutral 0.9 before collapsing conditions; known readouts require finite own-control CPM >5. Test this correction alone, a soft KD/0.8 reliability factor, and a >=50% gate for known measurements. Unknown readouts are not evidence of failed KD and are not hard-filtered. Original development covers unknown-KD targets CHKB, DHX35 and ZFP62; reserved confirmation includes additional such targets. Audit tables preserve every source row.
4. Keep V11's noise65 integer population emitter, own-target knockdown, expression gates and dual mean/bulk calibration unchanged. No combined response/emitter tuning is performed.

The first program screen used `median` instead of `nanmedian` for masked uncertainty. It was stopped after one K562 fold; those scores were not inspected for selection. The initial plan, exact initial code snapshot, and fold output are retained alongside the correction record. The repaired implementation was frozen before its outcomes.

## Validation and selection

All three summary screens exclude the entire evaluated source family. Summary normal-q proxies are insufficient for promotion: only the top two alternatives in each stage advance to actual generated-count tests. The robust stage was added during development after the program screen, and the quality stage followed the assay/KD audit; this is adaptive research, not a pre-registered comparison of all model classes.

The count comparison reuses the frozen V12 public target sets, controls, and two seeds, allowing exact V11 baseline replay. These are reused development data, not fresh independent biology. Predictions use 400 generated cells and genuine unique treatment rows (up to 400) with disjoint template/reference controls. Three biological source families and five technical strata are balanced equally within the hierarchy. Full independent references are primary. The sampled400 secondary panels retain and report missing NMAE rather than substituting zero.

The V13 policy requires complete primary panels, expression-error improvement, aggregate fidelity/reach drops no larger than 0.01, family drops no larger than 0.03, NMAE worsening no larger than 0.01, and PDS drop no larger than 0.02. Secondary available metrics receive the same guards. The greatest primary expression improvement among guard-passing candidates determines the winner. Official public source anchors are diagnostics, with the exact same truth panels and correct PDS target axis; small anchor denominators can amplify noise.

After development, freeze one candidate and evaluate it once with V11 on reserved count-confirmation targets. Do not retune after confirmation. If no candidate passes both stages, a user-requested new artifact may be packaged only as **exploratory**, with failed guards and tradeoffs explicit. Such an artifact is not a demonstrated improvement. Confirmation source-summary biology was previously inspected; only these generated-count outcomes were unopened.

The count scorer is the unchanged V12 CPU implementation verified against official cell-eval2 commit `0ce25cf4c495e76e60dfca4ddc8f822507aa0e18`. It uses raw counts, exact Scanpy normalization/log-rank arithmetic, tie-corrected no-continuity Wilcoxon, BH correction, reference-control CPM >5, the official NMAE eligibility gate, common PDS origin, and jackknife-corrected expression profiles. Private challenge anchors and treated ABC labels remain unavailable.

## Reproduce from project root (Windows Command Prompt)

```cmd
.\.competition_env\Scripts\python.exe collect_v13_rpe1.py
.\.competition_env\Scripts\python.exe audit_v13_rpe1.py
.\.competition_env\Scripts\python.exe screen_v13_responses.py
.\.competition_env\Scripts\python.exe test_v13_submission.py
.\.competition_env\Scripts\python.exe prepare_v13_validation.py
.\.competition_env\Scripts\python.exe evaluate_v13_submission.py --stage development --workers 4
.\.competition_env\Scripts\python.exe assess_v13_submission.py --stage development
.\.competition_env\Scripts\python.exe compare_v13_source_anchors.py
.\.competition_env\Scripts\python.exe screen_v13_robust.py
.\.competition_env\Scripts\python.exe test_v13_robust.py
.\.competition_env\Scripts\python.exe run_v13_robust_validation.py
.\.competition_env\Scripts\python.exe screen_v13_quality.py
.\.competition_env\Scripts\python.exe test_v13_quality.py
.\.competition_env\Scripts\python.exe run_v13_quality_validation.py
```

The production scripts preserve existing plans/models/packages; use a fresh checkout or recover the hashed records when reproducing an already completed run. Caches and downloaded data are excluded from Git; source provenance/hashes and result tables are tracked. `audit_v13_baseline_replay.py` checks all baseline profiles, jackknifes, and DE results against V12. Selection, confirmation, all-source fitting, generation, official packaging and independent payload replay are documented in `V13_SUBMISSION.md` after completion.

The quality stage uses a predeclared neutral 0.9 assumption for unquantifiable on-target RNA. A later audit recovered the author's construct-weighted KD metadata for these K562 targets; `published_knockdown_audit.json` records its distribution. This was found after the quality plan froze, so it is recorded as a future refinement and was not silently inserted into the tested model.

The V11 count emitter is unchanged. Its previously measured neutral-control significance and dispersion artifacts therefore remain a limitation; passing file checks or beating V11 locally does not establish null-calibrated biology.

Generation can reuse verified V11 count blocks only when the full response vector, original controls, population setting, and seed are exactly unchanged. The original H5AD is fully SHA-256 checked before reuse. Changed groups are regenerated. Preflight always generates fresh counts; the independent audit replays ABCD1, CHKB, DHX35, and ZFP62 in all three contexts, covering the missing-KD correction. This optimization preserves counts and does not create a new model variant.

Development selected `unknown_kd`: primary expression error decreases from 1.404989 to 1.402759 (about 0.16%), with small DE metric tradeoffs. Its first-seed public anchored overall decreases from 0.013840 to 0.013240, so the raw-guard pass is not an overall-score improvement. The Student-t alternative improves five primary metrics but worsens expression error by 0.002999 and fails that predeclared guard. Preserve both findings when interpreting the eventual artifact.

Reserved confirmation contains 30 target-seed cases from HCT116 and HEK293T only; H1 contributes 12 development cases and no confirmation cases. The full development set has 58 cases (12 H1, 24 HCT116, 22 HEK293T). Thus confirmation is a two-family check, not a three-family replication. It improves expression ratio by 0.035220 and fidelity by 0.001033, but fails the full-primary completeness guard because some panels have no eligible NMAE contrasts. No metrics or selection were changed to overcome this limitation.

A post-freeze arithmetic audit (`audit_v13_precision.py`) found 65 targets with corrected source KD and 235 others. The new pooling code retains float64 intermediates where V11 casts the pooled mean/variance to float32. For the other targets, the largest response difference is 8.345e-7 log2 units. This is included in the actual frozen/count-tested candidate, not corrected after confirmation. Overall 888/900 response vectors differ bitwise, so almost all groups must be regenerated. V13 is therefore a KD correction with small incidental precision changes, not a bitwise single-change ablation.
