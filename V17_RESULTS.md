# V17 learned response variation: measured results

**V17 is not released. V16 remains the recommended submission.**

We fitted target-specific response-strength variation from public H1/HCT116/HEK293T raw counts, excluding the entire evaluation family. The production-source fit supports nonzero excess variation for 85/300 targets. Five frozen arms were compared over 14 target/family pairs, two seeds and 400 cells per prediction: **140 generated populations**. Six additional real-control comparisons verified exact V16 null identity.

## Full-count development results

| Candidate | PDS | Fidelity | Reach | Jaccard | NMAE | Response-projection Wasserstein |
|---|---:|---:|---:|---:|---:|---:|
| v16 | 0.777778 | 0.606476 | 0.503430 | 0.022365 | 0.849526 | 1.215842 |
| learned_half | 0.777778 | 0.607336 | 0.482728 | 0.022403 | 0.859531 | 1.217988 |
| learned_full | 0.777778 | 0.606906 | 0.516763 | 0.022424 | 0.858893 | 1.215319 |
| generic | 0.777778 | 0.604744 | 0.422118 | 0.022229 | 0.863088 | 1.214686 |
| shuffled | 0.777778 | 0.607375 | 0.466070 | 0.022448 | 0.849286 | 1.213693 |

Higher PDS/fidelity/reach/Jaccard is better; lower NMAE and Wasserstein distance is better. These are family-balanced **public raw holdout metrics**, not official anchored scores. HCT116 NMAE is unavailable at both seeds because too few true DE genes pass the eligibility gate; NMAE averages the available H1/HEK293T panels only. No missing value is replaced by zero. Tiny control-replicate expression denominators make ratios unstable; the expression guard uses paired relative numerator changes.

## Decision and interpretation

- `learned_half`: rejected; failed checks: reach, nmae, family_reach, response_distribution_improves, beats_shuffled.
- `learned_full`: rejected; failed checks: beats_shuffled.

Full variation changes reach by +0.013333, fidelity by +0.000431 and Jaccard by +0.000059. NMAE worsens by +0.009367. Its relative expression-error improvement is only 0.00575%, and response-projection Wasserstein improvement is only 0.0430%.

The full candidate passes the available raw-score guards, but does not beat the target-shuffled variation control on two DE metrics. The improvement cannot be attributed reliably to learning the correct target-specific heterogeneity. Half-strength variation additionally fails reach, NMAE, worst-family reach and distribution checks. We retain the prospective policy rather than choose a favourable metric after seeing outcomes.

Reserved verification remains unopened when development fails. No new upload package is generated, and no leaderboard improvement is claimed. Existing V16 files and all numerical dependencies are preserved. The prepared V17 production adapters are gated and have not been exercised end-to-end without a supported finalist.

## Verification and limitations

- Generator invariants passed: exact V16 counts at zero response/zero variation; deterministic repeat; exact library sizes; nonnegative integer counts; moment constraints and bounded strength.
- Independent audit verified source-family exclusions, all actual/template/reference identities, pooled model arrays, both negative controls, cache hashes and metric-file hashes.
- Real-control null generation is bitwise V16; its average 1,946.5 generated DE calls versus zero actual control calls remains a baseline diagnostic problem. Adding response-strength heterogeneity does not repair that problem.
- Observational cell-state dependence and treatment-specific covariance are not identified by a scalar responder model. Excess source variance can reflect efficiency, survival, residual batch differences and other state changes.
- See [V17_STRATEGY.md](V17_STRATEGY.md) for literature, newly collected public metadata, input amendments and reproduction.

## What the new data discovery adds

The strongest new broad target source is the [2026 KOLF atlas](https://www.nature.com/articles/s41587-026-03199-w): 273 of our targets appear in its public browser. Browser mean NTC z-scores cannot establish cell-level covariance or be substituted for counts. Selective extraction of QC-filtered raw AnnData with matched controls is the justified next data step. Strong-perturbation-only filtering would bias response distributions and should not be used as the sole training population.

The [multi-donor iPSC release](https://doi.org/10.6084/m9.figshare.27989294) supplies additional biological diversity but sparse per-donor target counts. The targeted experiment contains five overlapping targets with 1,133–1,923 assigned cells pooled across donors, and only 96–144 in any single donor. Those are potential donor-held-out validation/prior-learning data, not 400-cell replicas of one context.

## Reproduction

```bat
.\.competition_env\Scripts\python.exe run_v17_experiment.py development --workers 6
.\.competition_env\Scripts\python.exe release_v17_submission.py
.\.competition_env\Scripts\python.exe plot_v17_evidence.py
.\.competition_env\Scripts\python.exe document_v17_submission.py
```

A finished experiment resumes from hash-verified caches. Collection/fitting plans refuse overwriting; see the saved plan and correction history for exact executed inputs.
