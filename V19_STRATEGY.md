# V19: covariance-conditioned response forces

Status: complete; failed the frozen release rule. No V19 package was generated. See `V19_RESULTS.md`. Prepared 2026-10-13 after PIE failed V18's frozen gates. This is a new algorithm from the post-V17 research plan, not a revised PIE confidence threshold.

## Hypothesis

Two contexts can have similar mean expression but different gene-program covariances and therefore different response susceptibilities. A response transferred as a fixed LFC vector does not represent this distinction. Approximate a source intervention as a force in a shared control-derived program space, then apply the query's control covariance. Query treatments are never observed.

This is inspired by [CIPHER](https://github.com/GoyalLab/CIPHER), whose inspected forward implementation fits an amplitude using measured query treatment shifts. V19 **does not use that oracle amplitude** and is not a reproduction of the published rank-one forward model.

## Model

For each held family, use only the other two raw-count families to fit a common basis. Source NTC rows are the same frozen left-control identities used by the existing experiments. Normalize per cell to CP10k using original full library totals, take natural log1p, subtract batch NTC means and regress log library size using only NTCs. Select the top 2,048 common-measured genes by mean adjusted control variance, excluding all 300 panel genes. Standardize with one training-derived scale and fit 16 shared SVD programs. Each source family contributes equal total weight to the basis.

Estimate source and query control covariance in this same program space. Use `0.8*C + 0.2*diag(C) + 0.05*mean(diag(C))*I`; do not independently whiten contexts to identity. Source intervention shifts are mean treatment program scores after subtracting their batch NTC expectation and NTC-estimated depth nuisance effect. Require at least 40 genuine source treatment cells and weight families by `n/(n+100)`.

For source `s` and query `q`, form `T = Sigma_q @ inverse(Sigma_s)`, bounding its singular values to [0.5, 2]. Apply `(T-I)` to the observed source program shift and average supported families. This is a context correction to V16's existing conserved response, not a replacement for all its source evidence.

Decode the correction through the shared program basis and training scale. Convert local mean-log-expression shifts to approximate log2 LFC using the query-control tangent `ln(2) * mean(CP10k/(1+CP10k))`, floored at 0.1. Clip each raw correction to +/-1, multiply by 0.25, apply the existing expression gate, and cap the combined LFC at +/-3. Own/panel coordinates never appear in the program basis; all unsupported genes retain V16. If fewer than 80% of basis genes are measured in the query, abstain.

This tangent is an approximation, not an exact inverse from mean log expression to mean CPM. Cross-assay covariance can include technical noise despite nuisance adjustment. The bounds and ablations are essential; source evidence does not establish that this model is causal.

## Prospective tests

Five arms: exact V16; covariance correction; diagonal source/query covariance; fixed query-program-label permutation; V16 with L1 amplitude matched to the covariance arm. All share the original V16 residual10 population generator, count seeds, own-target knockdown and 400-cell sampling. The development identities are explicitly reused from V18. They cover 14 target/family pairs, two seeds, and three full-family response exclusions.

Use the same available raw and worst-family guards. At least three metrics must improve including a DE metric. The real-covariance arm must beat **each** diagonal, shuffled and amplitude control on at least two of fidelity/reach/Jaccard. Missing NMAE stays missing. Only if these gates pass may the still-unopened eight corrected V17 target/family pairs be scored once under the same rule. A failed candidate is not retuned after verification.

Unit checks cover the identity/scalar covariance map, singular bounds, depth centering, zero-force fallback, correction bounds, missing-gene masking, low-coverage abstention and own-target preservation. Full count tests, exclusion/identity audits and full artifact validation are required before any release. The unchanged generator retains V16's null-DE and joint-variation limitations.

## Reproduction (Windows CMD)

```bat
.\.competition_env\Scripts\python.exe test_v19_covariance.py
.\.competition_env\Scripts\python.exe prepare_v19_experiment.py
.\.competition_env\Scripts\python.exe run_v19_experiment.py models
.\.competition_env\Scripts\python.exe audit_v19_model_inputs.py
.\.competition_env\Scripts\python.exe run_v19_experiment.py development --workers 6
.\.competition_env\Scripts\python.exe audit_v19_experiment.py development
.\.competition_env\Scripts\python.exe assess_v19_experiment.py development
```

Plan and fitted-model freezes refuse overwrites. Paired-case caches resume only with matching input hashes. Numerical arrays are ignored by Git; source hashes, plans, fit records, scored CSVs and audit records are retained. Consult `V19_RESULTS.md` for the executed result and release status.
