# Transfer validation audit

Completed 2026-10-07 using the frozen V7 public development profiles. This is a retrospective diagnostic of evaluation semantics, not a new model-selection result or a leaderboard estimate.

## Proven finding

The old local PDS subtracts a target-specific matched control profile from each prediction and truth profile. That turns identical predictions into target-specific effects. A label-independent prediction must instead score 0.5 under a common context control, as specified in the [official metric definition](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics.md).

All 94 saved NPZ profile hashes passed. The audit reproduced **all 162 historical model/panel PDS rows** to absolute tolerance 1e-12, establishing that the discrepancy comes from the control origin rather than an unrelated implementation change. Constant control-derived predictions score **0.854545-1.000000** historically across splits and seeds, versus **0.5** with the required common origin.

The table below averages the two stored seeds on the complete public panels. Historical and common values use identical saved predictions, truths, target lists, and scored gene masks. Only the PDS control origin changes.

| Held source | Targets | Constant: historical PDS | Constant: common PDS | V4 transfer: historical PDS | V4 transfer: common PDS |
|---|---:|---:|---:|---:|---:|
| H1 | 11 | 0.868182 | 0.500000 | 0.927273 | 0.668182 |
| HCT116 | 18 | 1.000000 | 0.500000 | 1.000000 | 0.612745 |
| HEK293T | 18 | 0.998366 | 0.500000 | 0.993464 | 0.553922 |

Near-perfect historical discrimination is therefore not evidence of near-perfect intervention prediction. The corrected origin does not prove that any previously rejected candidate should be adopted: expanded candidates also failed expression guards.

## Independent control-only comparison

`template_only` uses the actual 400 template controls for each saved test, with no downstream response and no synthetic count changes. The common expression ratio below uses their genuine delete-one jackknife, just as for the saved model predictions. Lower is better. Ratios are panel sums, with means across the two seeds; these are not scaled competition scores.

| Held source | Template-only common PDS | V4 common PDS | Template-only expression ratio | V4 expression ratio |
|---|---:|---:|---:|---:|
| H1 | 0.509091 | 0.668182 | 1.109541 | 1.058100 |
| HCT116 | 0.566993 | 0.612745 | 3.506195 | 4.057493 |
| HEK293T | 0.522876 | 0.553922 | 1.186216 | 1.471402 |

V4's expression ratio is worse than the template-only comparison in HCT116/HEK293T and modestly better in H1 under this diagnostic. This is consistent with difficulty transferring response magnitude, but the technical batch limitation below prevents treating these values as unbiased estimates of challenge-context biology. HCT116's weak corrected control-distance denominator also magnifies ratios.

## Exact procedure and safeguards

1. Verify the source-axis cache against its provenance hash; read its exact genes and targets.
2. Load the H1/HCT116/HEK293T selected raw counts with the existing loader, which verifies extraction file hashes and axes.
3. Reconstruct the original per-batch random split of non-targeting cells using seed 20261008. The template and reference populations are disjoint. Validate every saved template/reference cell ID against those populations.
4. Form one reference profile from the entire right-hand population in each source context: summed integer counts, normalized to 50,000, then log1p. Compute its delete-one jackknife with chunked, mergeable central moments to avoid materializing every dense cell at once.
5. Verify each saved profile hash. Evaluate all nine historical candidates, identical-profile controls, and genuine template-only controls. Exclude all 300 official target coordinates for PDS.
6. Reproduce old PDS rows; recompute common-control PDS. Recompute expression with a common denominator origin and whole-gene scalar jackknife corrections. Retain the signed numerator and report unresolved denominators explicitly instead of flooring them to 1e-12. No denominator in this run was nonpositive.
7. Record profile/sidecar hashes, historical dependencies, current numerical dependency hashes, reference/template cell IDs, source records, script hash and output hashes in the audit manifest.

The control populations contain 768/768 left/right H1 cells across 48 batches, 1,742/1,742 HCT116 cells across 109 batches, and 3,558/3,558 HEK293T cells across 223 batches. This audit does not expand these selected-cell caches into full-release training data.

Self-tests passed: identical-prediction PDS 0.5, zero-effect PDS 0.5, perfect-profile PDS 1, joint-label permutation invariance, rejection of target-specific control arrays, chunked-versus-original jackknife equivalence, zero-sampling expression limits, and unresolved-denominator reporting. A synthetic counterexample reproduces historical PDS 1 with identical predictions.

## What this audit cannot establish

- The public panels have 11/18/18 targets, not the competition's 300 targets, and lack official anchors.
- A single pooled public control may mix technical batches. The cached predictions were made from templates matched to the actual treated batch proportions; that input is unavailable in the challenge. `template_only` reveals this nuisance baseline but does not remove it.
- Calibration and confirmation targets were already inspected during V7. This retrospective audit must not be used to select a winner and describe it as independently confirmed.
- DE metrics require per-cell counts and a shared control gene universe. They are not recomputed from saved profile arrays; prior target-matched DE metrics are also not exact competition-context replicas.
- No response is refitted, no V7 candidate is reselected, and no published payload or frozen numerical module is altered. The corrected PDS function is available in the new root audit script; historical `v6_methods.py` remains preserved for reproduction.

## Reproduction and outputs

From Windows Command Prompt:

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe audit_transfer_validation.py --self-test
.\.competition_env\Scripts\python.exe audit_transfer_validation.py
```

The full audit requires the existing ignored V7 profile NPZs and selected public raw-count caches. Reconstruct those with the pinned procedures in [V7_MODEL_STRATEGY.md](V7_MODEL_STRATEGY.md) and [V7_SUBMISSION.md](V7_SUBMISSION.md) when they are absent.

Outputs: [profile_audit.csv](analysis/transfer_validation_audit/profile_audit.csv), [profile_summary.csv](analysis/transfer_validation_audit/profile_summary.csv), [audit_manifest.json](analysis/transfer_validation_audit/audit_manifest.json). Research interpretation and the replacement validation design are in [VCC2026_RESEARCH_AND_GAPS.md](VCC2026_RESEARCH_AND_GAPS.md).
