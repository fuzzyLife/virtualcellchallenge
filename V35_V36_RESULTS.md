# V35–V36 results and submission

V24 remains the published reference at overall 0.14299785521745675. These are reused public development tests, not private anchored leaderboard scores. A new experimental package was explicitly requested and prepared; it is not a confirmed better model.

## Gene-level summary response tests

| Model | Squared LFC error change vs V24 | Absolute LFC error change | Strict summary gate |
|---|---:|---:|---|
| V35 | -0.192% | -0.323% | Fail |
| V36 | +0.014% | +0.125% | Fail |

Whole H1, whole Orion (HCT116+HEK293T), and whole CD4 exclusions are averaged equally. The learned retrieval is compared with response identities shuffled within each source and fixed attention at inference. These controls are ablations, not separately trained negative-control models. V35 learns detailed experimental retrieval and a bounded residual; V36 adds independent raw-population views and a dispersion head. Neither is an exact reproduction of the cited papers.

## Coverage

| Excluded family | Shared modeled output genes | Eligible raw population groups, fitting and validation |
|---|---:|---:|
| H1 | 5595 | 616 |
| Orion | 5540 | 73 |
| CD4 | 6221 | 593 |

Population groups repeat interventions across contexts; repeated sampled views are not independent biological experiments. In particular, the small Orion-excluded population branch limits evidence for broad transfer. Other genes retain V24 transfer-effect values exactly; normalization and sampling can change their emitted counts.

## Actual-cell mean tests

| Model | Squared LFC error change | Absolute LFC error change | CPM mean squared error change |
|---|---:|---:|---:|
| V35 | -0.551% | -0.941% | -5.531% |
| V36 | +0.097% | -0.543% | -4.664% |
| TARGET_CENTERED25 | -0.614% | -0.555% | -2.268% |
| RETRIEVAL_CONTRAST25 | -0.058% | -0.057% | -0.176% |

All 44 previously spent raw mean cases were evaluated, with the matching family's responses excluded. Different denominators mean this table and the summary table are not interchangeable. No fresh held-out lineage or D3 donor result is claimed.

## Raw count diagnostics

| Arm | PDS change | FID change | REACH change | JAC change | NMAE change | Relative expression error change |
|---|---:|---:|---:|---:|---:|---:|
| v35 | -0.180556 | -0.010722 | -0.052271 | -0.001940 | +0.014509 | +5.739% |
| v36 | -0.208333 | -0.011052 | -0.043528 | -0.001971 | +0.020037 | +18.435% |
| v36_dispersion | -0.208333 | -0.010323 | -0.054874 | -0.002041 | +0.019112 | +18.435% |
| v35_population_dispersion | -0.180556 | -0.010653 | -0.039555 | -0.002318 | +0.018911 | +5.739% |
| target_centered25 | -0.041667 | -0.002304 | +0.001181 | -0.001769 | +0.002895 | -2.537% |
| retrieval_contrast25 | +0.000000 | -0.003645 | +0.008023 | -0.000661 | +0.001437 | -0.515% |

PDS/FID/REACH/JAC increases are favorable; NMAE/expression-error decreases are favorable. Twenty treatment cases: first four cases from each of five fixed public contexts. D1 treatment cases use 200 cells and are lower-power diagnostics; other cases use 400. V35/V36 means first use unchanged V24 count emission; the final arm separately tests learned gene-wise dispersion weights. Undefined values are never replaced with zero. Metrics average defined contexts within each study; unresolved studies fail their guard.

| Context | Actual null DE calls | V24 | V35 | V36 | V36 dispersion |
|---|---:|---:|---:|---:|---:|
| CD4_Rest_D1 | 1 | 516 | 516 | 516 | 516 |
| CD4_Stim48_D2 | 0 | 753 | 753 | 753 | 753 |
| H1 | 0 | 1415 | 1415 | 1415 | 1415 |
| HCT116 | 1 | 780 | 780 | 780 | 780 |
| HEK293T | 0 | 355 | 355 | 355 | 355 |

Null fidelity is not treatment recovery. Existing null artifacts are explicitly retained in the report; a plausible-looking cell cloud is not sufficient validation.

## Selected artifact

- Selected arm: **retrieval_contrast25**.
- Status: **exploratory candidate; strict improvement guards not all passed**.
- Selected mean error changes versus V24: squared LFC -0.058%, absolute LFC -0.057%, CPM -0.176%.
- Strict public improvement guards: **False**.
- New blind confirmation: **False**. D3 reserve opened: **False**.
- Package: `experiments/v36_population_views/prediction.vcc`.
- Bytes: 4,073,666,560; SHA-256: `f9e71ec9f4608603139acf6ecfac03025ef3e0c76b0c0a6eec503f9011c04dd7`.
- Official prep/validation: passed. Automatically uploaded: False.
- The original V24 package and champion record remain unchanged.

The secondary V35-mean/V36-dispersion component check was added after summary, raw-mean and partial count outcomes were visible. It uses identical cases, fixed weights and a separately frozen rule, and remains explicitly exploratory. Original primary selection is preserved in `primary_selection.json`.

A second exploratory correction was added after the full primary outcomes exposed a large PDS/expression regression. Two fixed quarter-strength transforms were tested: remove the target-shared correction, or retain genuine-minus-shuffled retrieval dependence. Complete 300 predictions define centering; no treatment labels enter the transforms. The same 44 means and 20 count cases were rerun for each. `precontrast_selection.json` preserves the prior decision; contrast results and audits remain separate.

The pre-outcome primary rule chooses between the two experimental mean candidates using equal-study raw mean error, preferring count noninferiority when available. Dispersion is adopted only if its separate predefined comparison passes. If strict gains are absent, the package remains exploratory; no leaderboard improvement is asserted.

## Verification and reproduction

Seven boundary tests passed. All 6 outer models' source/family/target exclusions and saved summary scores were audited. A separately written NumPy implementation checks the neural forward. Explicit disjoint cell identities are checked for each source/target's first pair; all pairs are asserted during construction. All 100 primary count-arm DE-set summaries and 20 primary treatment PDS/expression panels were replayed independently. The exploratory combination and contrast arms add 60 DE-set summaries and 15 profile panels, also independently replayed from saved arrays. Reach's p-value tie order and complete raw Wilcoxon computation were inherited from the frozen scorer, not independently regenerated.

Production preflight caught and repaired object-string label serialization before count generation. Initial files and generator source are preserved; numeric predictions were verified exactly unchanged without unpickling. The corrected adapter and independent production effect audit are recorded. Full count/library/gene/order/density validation and official lossless packaging are recorded in analysis receipts. A separate final audit independently checks every one of the 360,000 row libraries and its sampling-provenance/obs alignment. Large model/data/prediction arrays are ignored in Git; root scripts and metadata document their regeneration. See [strategy and exact commands](V35_V36_STRATEGY.md), [safe numerical replay](reproduce_v35_v36.py) and `requirements-v35-science.txt`.

Original plans and code are preserved. A separate pre-outcome selector amendment corrected the scorer field name from a nonexistent expression-error label to `common_expression_ratio`, without changing thresholds or the comparison. No numeric result was available when that initial field-name correction was made. A later separate contrast-selector repair handles undefined HCT116 NMAE using the unchanged primary aggregation; original predictions, failed selector and logs are preserved. One fixed training seed and one fixed emission seed per case were used; small gains are not statistical proof of generalization.

Windows Command Prompt upload command:

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe -m vcc submit .\experiments\v36_population_views\prediction.vcc -m "Retrieval-contrast response correction v36" --wait
```
