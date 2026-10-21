# V89 measured results

**Rejected; no submission package generated.**

The candidate is fixed by the prospectively recorded protocol. All public biology is reused, both technical seeds share donor biology, and no hidden score improvement is established. V66/V82 remains the best published reference until a new hosted result exists.

## Broad public means

655 target/context pairs with equal study votes. Lower profile error is favorable; higher raw PDS is favorable. These values are not anchored competition scores.

| Arm | Profile error / V66 | Error reduction | PDS change |
|---|---:|---:|---:|
| interval25 | 0.996265839 | 0.3734% | -0.000536222 |
| shared25 | 0.995987093 | 0.4013% | -0.000752674 |
| interval25_shuffled | 0.993357747 | 0.6642% | -0.009481510 |

Strong broad checks: profile=False, worst_study=True, pds=True, beats_shared=False, beats_shuffle=False.

## Paired integer-count results

| Seed | Arm | Profile MSE | PDS | Fidelity | Reach | Jaccard | NMAE |
|---|---|---:|---:|---:|---:|---:|---:|
| 0 | v66 | 0.006497002 | 0.565934066 | 0.518864975 | 0.006449255 | 0.152418328 | 0.983722336 |
| 0 | interval25 | 0.006478115 | 0.565934066 | 0.521369712 | 0.008381768 | 0.153310556 | 0.981878256 |
| 0 | shared25 | 0.006478978 | 0.565934066 | 0.516549449 | 0.006998598 | 0.152317605 | 0.983218426 |
| 0 | interval25_shuffled | 0.006467136 | 0.560439560 | 0.518889108 | 0.007353236 | 0.152470531 | 0.983920272 |
| 620000 | v66 | 0.006444256 | 0.549450549 | 0.535880857 | 0.007216457 | 0.146324904 | 0.983286876 |
| 620000 | interval25 | 0.006424781 | 0.549450549 | 0.536206105 | 0.007216457 | 0.144698378 | 0.982504455 |
| 620000 | shared25 | 0.006425827 | 0.549450549 | 0.539098924 | 0.006839837 | 0.145362871 | 0.982698121 |
| 620000 | interval25_shuffled | 0.006414281 | 0.549450549 | 0.534846078 | 0.007535618 | 0.145267808 | 0.983086284 |

## Frozen decision

- Strong support: False.
- Bounded experimental release allowed: False.
- Broad safety checks: profile=True, worst_study=True, pds=True.
- Seed 0 experimental checks: profile_MSE=True, pds=True, fidelity=True, jaccard=True, nmae=True, two_gains=True.
- Seed 0 improved endpoints: profile_MSE, fidelity, reach, jaccard, nmae.
- Seed 0 strict checks: profile_MSE=True, pds=True, fidelity=True, reach=True, jaccard=True, nmae=True, three_gains=True, beats_shared25=True, beats_interval25_shuffled=False.
- Seed 620000 experimental checks: profile_MSE=True, pds=True, fidelity=True, jaccard=False, nmae=True, two_gains=True.
- Seed 620000 improved endpoints: profile_MSE, fidelity, nmae.
- Seed 620000 strict checks: profile_MSE=True, pds=True, fidelity=True, reach=True, jaccard=False, nmae=True, three_gains=True, beats_shared25=True, beats_interval25_shuffled=False.

V89 releases the earlier V86 interval25 configuration; it is not a newly learned architecture. Public predictions and all56 first-seed records were adopted from V86, after exact hash checks. The second-seed test generated42 new populations and reused14 V84 baselines. The known first-seed donor result informed candidate selection, so the repeat test is technical verification on reused biology, not a blind clinical or biological confirmation.

The stronger criteria still matter: the candidate did not establish superiority to all controls. The bounded experimental rule does not establish a private score gain. The V86 finite-iteration fit and normal-moment interval remain approximations; production source-model metadata reports actual iteration counts without claiming convergence.

Root code and compact evidence are versioned; large arrays, models, cells and packages remain local and ignored. See V89_STRATEGY.md and the frozen JSON protocols for the full reproduction sequence.
