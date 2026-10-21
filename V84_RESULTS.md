# V84 results: measured cross-context response operators

**New submission qualification: False.** Development selected `None`. Both fixed operators failed the original expression-improvement requirement. The primary was nevertheless evaluated on the separate donor as prespecified. No new upload package was generated; V66/V82 remains the best published project result, overall **0.1508821875498449**. No new hosted score is claimed.

[Frozen strategy and reproduction](V84_STRATEGY.md) | [Research and project gaps](POST_V82_RESEARCH.md) | [Coverage correction](V83_COVERAGE_AUDIT.md)

## Source-task results need strong baselines

Train on 513 shared essential interventions; evaluate 125 held interventions on 6,087 nonpanel outputs. No competition target is in that training panel. Lower ratios are better. Copying K562 is the identity baseline; the constant baseline is the destination fitting-set average response.

| Operator | Destination | Error / copied K562 | Error / constant response |
| --- | --- | ---: | ---: |
| program25 | RPE1 | 0.718888 | 1.125306 |
| program25 | hepg2 | 0.793689 | 1.107489 |
| program25 | jurkat | 0.913422 | 1.385339 |
| diagonal25 | RPE1 | 0.803621 | 1.257941 |
| diagonal25 | hepg2 | 0.757951 | 1.057621 |
| diagonal25 | jurkat | 0.658055 | 0.998038 |

The primary improves on copying K562 but loses to the constant-response comparator in all three destinations. The diagonal operator only narrowly improves on the constant comparator in Jurkat. These results do not establish recovery of transferable perturbation identity. The operators were not refitted on held targets.

## Reused public development screen

All 655 target/context pairs were included, with legacy recipient-study response evidence excluded. Macro averages give H1, Orion and CD4 equal study weight. These are reused public cases, not independent confirmation.

| Candidate | Profile error / V66 | PDS change | Pass every original development criterion |
| --- | ---: | ---: | --- |
| program25 | 1.000047961 | 0.002195 | False |
| diagonal25 | 1.000761708 | -0.014450 | False |
| program25_shuffled | 1.000524580 | -0.000752 | negative control |
| diagonal25_shuffled | 1.000397765 | -0.013800 | negative control |

Require profile ratio <=0.995, worst-study ratio <=1.03, PDS loss <=0.005, and better profile error than the separately fitted shuffled model. No arm qualifies. The primary passes the latter three requirements but leaves expression error effectively unchanged. The diagonal arm also fails discrimination and the shuffled comparison.

## Project-unseen D3 donor: final emitted counts

The 14 targets were selected from metadata before scoring. The whole CD4 response family and merged donor teachers were excluded. Query inference used 5,006 NTCs; all comparisons used a disjoint 4,994-cell reference and matched 400-cell actual/template populations. Each technical seed is assessed separately. The donor belongs to a previously used study; the seeds may overlap cells and are not independent biological replicates.

Raw public metrics are shown below. Smaller profile error and NMAE are better; larger PDS, fidelity, reach and Jaccard are better. These are not hosted anchor-scaled scores or a predicted overall score.

| Seed offset | Arm | Profile error | PDS | Fidelity | Reach | Jaccard | NMAE |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | v66 | 0.006497 | 0.565934 | 0.518865 | 0.006449 | 0.152418 | 0.983722 |
| 0 | program25 | 0.006503 | 0.549451 | 0.521334 | 0.005214 | 0.150995 | 0.983727 |
| 0 | program25_shuffled | 0.006501 | 0.571429 | 0.518168 | 0.007007 | 0.149657 | 0.981962 |
| 620000 | v66 | 0.006444 | 0.549451 | 0.535881 | 0.007216 | 0.146325 | 0.983287 |
| 620000 | program25 | 0.006447 | 0.543956 | 0.534075 | 0.007491 | 0.146124 | 0.983426 |
| 620000 | program25_shuffled | 0.006444 | 0.565934 | 0.534987 | 0.006632 | 0.145348 | 0.986486 |

| Seed | Primary profile error change | Primary PDS change | Failed original donor checks |
| ---: | ---: | ---: | --- |
| 0 | +0.0895% | -0.016484 | profile, pds, jaccard, beats_shuffle, three_gains |
| 620000 | +0.0436% | -0.005495 | profile, pds, beats_shuffle, three_gains |

The comparator is the V66 algorithm recomputed with all CD4 responses excluded, not the submitted ABC array. Complete per-target records, NMAE availability lists and individual guard decisions remain in the two donor assessment JSON files. No target, strength or threshold was selected from donor outcomes.

## Conditional diagnostic: mean versus population shape

After the fixed primary failed, the presaved protocol selected ADNP, MAP2K2, MYRF, PCGF1, USP28 and ZNF25. Each target uses disjoint 200-cell treatment fitting/scoring halves. Treatment labels enter only these non-deployable diagnostic populations.

| Diagnostic population | Profile error | Fidelity | Reach | Jaccard | NMAE | CPM contract holds for all six |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| baseline | 0.009374 | 0.498093 | 0.008803 | 0.107467 | 0.983400 | True |
| fitted_mean_old_shape | 0.006119 | 0.828716 | 1.000000 | 0.160845 | 0.406014 | True |
| predicted_mean_fitted_shape | 0.009374 | 0.103482 | 0.005669 | 0.047619 | 0.984267 | False |
| real_fit | 0.006119 | 0.929664 | 1.000000 | 0.314528 | 0.406024 | True |

Within each mean pair, gene totals, row totals and raw profiles match exactly: `baseline` versus `predicted_mean_fitted_shape`, and `real_fit` versus `fitted_mean_old_shape`. The independent audit replays profiles, BH and direction metrics across all 24 populations and checks 192 original tied-rank coordinates. A failed CPM contract means that a projected population cannot be treated as isolating shape on every statistic.

Across mean pairs, fitted treatment depth and projection behavior also differ. Measured fitting-half means are noisy; six targets and 200-cell scoring halves do not establish a universal bottleneck. Oracle populations cannot train or select an upload model.

## What this changes in our priorities

On this limited diagnostic, supplying measured fitting-half means reduces profile error by 34.72% and raises raw fidelity from 0.498093 to 0.828716. That is evidence of response-mean headroom in these six cases, not evidence that our available inference inputs can recover those means.

DE overlap still rises from 0.160845 with fitted means and the old shape to 0.314528 using real fitting-half cells. Matching means therefore leaves distributional errors. The fitted-shape-at-baseline-means branch fails its CPM contract, so its poor scores cannot isolate a pure shape effect.

The next justified modeling target is accurate target-specific expression change in the scoring statistics, demonstrated on a whole held-out study and against constant/shuffled controls before adding a new population generator. More source data must first pass explicit target coverage or paired-response transfer tests. Another blend or dispersion setting without that evidence is unlikely to resolve the observed gap.

## Integrity and preserved failures

- Both donor seeds have 42 comparison populations, 84 total, with 672 independently checked tied-rank coordinates. Independent replay verifies BH, profiles, panel PDS and DE direction metrics.
- The exclusion audit checks 28 lane-matched views and bitwise invariance after poisoning all excluded CD4 responses, errors, doses, masks and source profiles. Access to the CD4 teacher is forbidden in the replay.
- Independent operator checks verify the normal equations, fitting split, scales, diagonal formula, held-target errors, zero-input identity and frozen source/code hashes.
- V83 direct augmentation stopped at zero target coverage. The earlier V17 finding and initial parser engineering corrections are retained. Missing targets were never manufactured as biological zero votes.
- The scDEBART code audit checks reviewed attention classes and configuration at a pinned commit. No released checkpoint was executed, and no claim is made about the code used for its published benchmarks.
- Large source files, arrays, count populations and logs are ignored by Git. Compact provenance, code, fixed plans, per-case scores and audits are tracked. Original scientific releases were not rewritten.

## Reproduce and interpret

Use the commands in `V84_STRATEGY.md`, then run `report_v84_results.py`. The execution protocol records the model/code/input hashes before recipient evaluation. `finish_v84_donor.py` resumes completed case files only after checking their array hashes.

The new paired-evidence approach was implemented and tested rather than assumed to work from a paper. It does not currently justify replacing V66. See `POST_V82_RESEARCH.md` for the evidence gap and the distinction between our task and published benchmarks.
