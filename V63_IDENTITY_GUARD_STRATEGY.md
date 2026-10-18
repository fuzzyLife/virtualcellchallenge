# V63: protect intervention identity during assay-space calibration

V62's new assay-space response representation improves raw mean error across all three held studies, but its 60 count comparisons fail the fixed PDS and Jaccard limits. Those negative results remain in V62_RESULTS.md. V63 is a separately frozen, outcome-informed follow-up; it is not a claim that the guard was proposed before seeing those failures.

## New hypothesis

The V24 model already carries useful target identity. A correction that improves expression error can rotate that signal and lose discrimination. V63 accepts the fixed V62 correction only as far as it preserves the original model's distinctions between targets. It does not fit another neural network or use the held recipient's treated cells.

For each query, construct the original V24 pooled expression response for all 300 targets. Remove all 300 panel genes from the geometry. Let d0[t] be the response and u[t] its unit vector. The retained margin against any other target s must satisfy:

```
dot(new_response[t], u[t] - u[s]) >= 0.8 * dot(d0[t], u[t] - u[s])
```

The response refers to log1p(50000 * pooled gene proportion) minus the query-control profile. The fixed estimator tries the complete V62 correction, then halves its interpolation fraction until every constraint holds. It interpolates normalized pooled probabilities, converts back to the exact own-knockdown convention, and checks the realized float32 effects after closure. Original V24 is always a feasible fallback. The 80% retention constant is frozen before V63 count scores; no threshold is tuned per held context.

This is a check against V24's predicted templates, not measured truth. Preserving old margins cannot guarantee preservation of true PDS. That remains a measured test. The method is an independently derived safeguard motivated by the failure observed in V62, building on the assay-space/source-common experiment and its cited response-decomposition research.

## Tests and boundaries

The source estimator, original templates, controls and all target labels come from permitted inputs. Complete H1, Orion and CD4 exclusions match V62. V62's adversarial held-study poisoning audit still applies to the unchanged source estimator. New synthetic tests verify probability closure, the own-target factor convention, and every target-margin constraint; production receives an independent per-target inequality audit.

One fixed candidate, guarded_assay, is tested against V24 and the equally guarded cyclic source-target-shuffled diagnostic. Development uses the same 20 public target populations. Identical V24 baseline profile/DE statistics are copied byte-for-byte from V62 with receipts, because source populations, seeds and generator are identical. Count logs retain the inherited evaluator's V62 print prefix; V63 result paths and receipts identify the actual experiment. This is computational reuse, not new independent evidence.

The original confirmation policy permits remaining cases only after development support. In this actual run development failed; run_v63_additional_counts.py instead labels the 24 remaining configurations at two seeds as post-failure public diagnostics. They do not provide confirmation credit or change the fixed candidate. These public studies and their mean outcomes were used before; this is case/seed confirmation, not a blind biological reserve. D3 remains unopened. The unchanged generator's excess null DE remains unresolved. The broader identity diagnostic includes all previously extracted targets with at least 30 real cells, with actual pooled-count truth; it is a desired-mean diagnostic before integer rounding and cannot be reported as a full 300-target leaderboard test.

## Reproduction

Use the existing Python 3.12 scientific environment, exact saved V7/V24/V34/V37 inputs and the V62 source estimator. Original plans refuse overwrite.

```cmd
SCIENCE_PYTHON test_v63_guard.py
SCIENCE_PYTHON v63_identity_guard.py freeze
SCIENCE_PYTHON v63_identity_guard.py prepare
SCIENCE_PYTHON run_v63_development_parallel.py
SCIENCE_PYTHON run_v63_additional_counts.py
SCIENCE_PYTHON audit_v64_unique_panels.py --version 63 --stage additional
SCIENCE_PYTHON analyze_v62_v63_diagnostics.py
SCIENCE_PYTHON diagnose_v63_amplitude.py
SCIENCE_PYTHON audit_v64_unique_mean_control.py
SCIENCE_PYTHON evaluate_v63_broad_identity.py
```

The parallel development driver waits for prediction receipts and uses at most two independent count workers. Alternatively run development separately per source and assess-development afterwards. Production scripts are fit_v63_submission.py, build_v63_submission.py, package_v63_submission.py, verify_v63_effects.py and validate_v63_submission.py. Fit and build only the locked candidate with a documented count_assessment.json. A failed candidate requires explicit experimental mode and must not be called a validated improvement. Only the user uploads; no Git push is performed.

Validation correction: the 44 configurations span 30 unique target/context combinations and 34 distinct measured row sets; repeated sets are not new biological evidence. Additional HEK293T scoring receives canonical target deduplication in a separate audit, preserving the original failed results. See V62_V64_VALIDATION_AUDIT.md.
