# V22 measured results

V22's bounded state-dependent response model is **rejected by its frozen development policy**. V16 remains the published champion. V23 is a separately frozen mean-preserving extension.

Completed45 oracle populations (nine disjoint200/200 calibration/truth contrasts),196 full-count development populations (14 public target/family pairs, two seeds, seven arms), six full-axis null populations and five held-family/all-data factor fits. Independent audit verifies source exclusions, cases, caches and all28 exact prior V16 baseline count hashes.

## Development

Raw metrics: higher PDS/fidelity/reach/Jaccard, lower NMAE/expression error/shape distance are better. HCT116 NMAE is missing. The expression ratio includes an unstable HCT116 denominator; release uses paired absolute numerators.

| candidate | pds | fidelity | reach | jaccard | nmae | common_expression_ratio | shape_wasserstein |
| --- | --- | --- | --- | --- | --- | --- | --- |
| amplitude | 0.777778 | 0.608342 | 0.498703 | 0.0221933 | 0.859107 | 67.9864 | 0.347698 |
| conditional | 0.772222 | 0.609671 | 0.507362 | 0.0226227 | 0.864304 | 67.7322 | 0.340442 |
| conditional_half | 0.777778 | 0.609319 | 0.452425 | 0.0225337 | 0.859391 | 67.7263 | 0.347323 |
| independent | 0.772222 | 0.60766 | 0.491339 | 0.0222404 | 0.858919 | 67.7174 | 0.355188 |
| scalar | 0.777778 | 0.608751 | 0.518298 | 0.0220843 | 0.856969 | 67.7218 | 0.345628 |
| shuffled_state | 0.772222 | 0.607404 | 0.440856 | 0.0223384 | 0.868869 | 67.7201 | 0.345526 |
| v16 | 0.777778 | 0.606476 | 0.50343 | 0.0223651 | 0.849526 | 67.724 | 0.348377 |

Conditional-minus-V16: fidelity+0.0031950, reach+0.0039321, Jaccard+0.0002577, shape distance-0.0079345 and paired expression error-0.0100%. It beats scalar, independent and shuffled-state DE controls. However, NMAE worsens+0.0147788 above the frozen+0.01 allowance. This failure rejects release; the threshold was not relaxed. PDS falls-0.0055556 within its guard.

## Factorial oracle

| candidate | pds | fidelity | reach | jaccard | nmae | common_expression_ratio | shape_wasserstein |
| --- | --- | --- | --- | --- | --- | --- | --- |
| measured_mean_compact | 0.777778 | 0.621517 | 1 | 0.0292038 | 0.278802 | 2.77432 | 0.306034 |
| measured_mean_treated_shape | 0.777778 | 0.92847 | 1 | 0.598701 | 0.278245 | 2.27921 | 0.21315 |
| pred_mean_compact | 0.722222 | 0.561466 | 0.682365 | 0.0138494 | 0.895663 | 1.93757 | 0.318013 |
| pred_mean_treated_shape | 0.722222 | 0.415346 | 0.721143 | 0.269888 | 0.890623 | 1.94176 | 0.312953 |
| real_replica | 0.777778 | 0.92847 | 1 | 0.598701 | 0.278248 | 2.27919 | 0.213568 |

Measured averages plus measured shape reach Jaccard0.5987 versus0.01385 for predicted averages/compact. Mean-only0.02920; shape-only0.26989 but fidelity worsens. These reuse nine public contrasts,200-cell treatment splits and common calibration-treatment depths. NMAE has limited eligibility. They are diagnostic ceilings, not deployable models or estimates of private scores. The measured-shape/mean arm nearly reproduces its calibration-cell replica by construction.

## Null and data audit

Zero-effect V22 matches V16 bitwise. Both still produce about1,487 DE calls versus zero calls from the real-null samples; identity establishes non-regression, not good null calibration. Guide efficacy and causal basal state remain incompletely identified.

Donor1 CD4 Rest was downloaded/MD5-verified, but only two competition targets reach400 cells. Metadata-only probes find donor2 Rest15, donor1 Stim48hr2 and donor4 Rest1. Donor2 was selected and independently MD5-verified. Its retained cells are all low_quality=False;23 lanes are matched. Original donor1 extraction is archived. The prospective source amendment records lane_id mapping and the powered-source choice. No V22 reserve scores were opened after rejection. V23 freezes the same eight metadata-selected donor2 targets before reserve scoring.

A subsequent no-DE diagnostic reduced V22-to-V16 mean-CPM L1 drift from0.00201614 to0.00000206694 using integer row/column-preserving swaps. This motivated V23, with unchanged V22 factors and unchanged rejection. Methods/commands: V22_STRATEGY.md. Evidence: experiments/v22_state_response/analysis/.
