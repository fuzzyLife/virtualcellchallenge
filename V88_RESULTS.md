# V88 measured results

**Rejected; no submission package generated.**

The candidate is fixed by the prospectively recorded protocol. All public biology is reused, both technical seeds share donor biology, and no hidden score improvement is established. V66/V82 remains the best published reference until a new hosted result exists.

## Broad public means

655 target/context pairs with equal study votes. Lower profile error is favorable; higher raw PDS is favorable. These values are not anchored competition scores.

| Arm | Profile error / V66 | Error reduction | PDS change |
|---|---:|---:|---:|
| exact125 | 0.995757940 | 0.4242% | -0.001833226 |
| shared125 | 0.996073075 | 0.3927% | -0.001303294 |
| exact125_shuffled | 0.994180175 | 0.5820% | -0.008356000 |

Strong broad checks: profile=False, worst_study=True, pds=True, beats_shared=True, beats_shuffle=False.

## Paired integer-count results

| Seed | Arm | Profile MSE | PDS | Fidelity | Reach | Jaccard | NMAE |
|---|---|---:|---:|---:|---:|---:|---:|
| 0 | v66 | 0.006497002 | 0.565934066 | 0.518864975 | 0.006449255 | 0.152418328 | 0.983722336 |
| 0 | exact125 | 0.006484047 | 0.560439560 | 0.520316603 | 0.006732848 | 0.151519872 | 0.980417142 |
| 0 | shared125 | 0.006487074 | 0.565934066 | 0.519441577 | 0.004008354 | 0.151011046 | 0.984596553 |
| 0 | exact125_shuffled | 0.006469921 | 0.571428571 | 0.519979129 | 0.007176655 | 0.152907499 | 0.981763486 |
| 620000 | v66 | 0.006444256 | 0.549450549 | 0.535880857 | 0.007216457 | 0.146324904 | 0.983286876 |
| 620000 | exact125 | 0.006429655 | 0.549450549 | 0.538072879 | 0.005791583 | 0.145438393 | 0.984132599 |
| 620000 | shared125 | 0.006434267 | 0.549450549 | 0.537118223 | 0.007135887 | 0.146197953 | 0.981608166 |
| 620000 | exact125_shuffled | 0.006415094 | 0.554945055 | 0.537175670 | 0.006313293 | 0.146625790 | 0.983757226 |

## Frozen decision

- Strong support: False.
- Bounded experimental release allowed: False.
- Broad safety checks: profile=True, worst_study=True, pds=True.
- Seed 0 experimental checks: profile_MSE=True, pds=False, fidelity=True, jaccard=False, nmae=True, two_gains=True.
- Seed 0 improved endpoints: profile_MSE, fidelity, reach, nmae.
- Seed 0 strict checks: profile_MSE=True, pds=False, fidelity=True, reach=True, jaccard=False, nmae=True, three_gains=True, beats_shared125=True, beats_exact125_shuffled=False.
- Seed 620000 experimental checks: profile_MSE=True, pds=True, fidelity=True, jaccard=False, nmae=True, two_gains=True.
- Seed 620000 improved endpoints: profile_MSE, fidelity.
- Seed 620000 strict checks: profile_MSE=True, pds=True, fidelity=True, reach=True, jaccard=False, nmae=True, three_gains=False, beats_shared125=True, beats_exact125_shuffled=False.

The converged fit and numerical mixture intervals did not meet the overlap limits on either seed; first-seed PDS also failed. No control or alternative strength was promoted. V88 produced84 new populations plus28 reused baselines across its two seeds, with672 fresh tied-rank coordinates and224 previously checked baseline coordinates. Independent source, profile and count-statistic audits passed. Four production adapters were written in advance but were not run because release was blocked.

Root code and compact evidence are versioned; large arrays, models, cells and packages remain local and ignored. See V88_STRATEGY.md and the frozen JSON protocols for the full reproduction sequence.
