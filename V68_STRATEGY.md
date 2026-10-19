# V68: predict populations by balancing and selecting real cells

V67's new intervention-network model gave only a 0.0018% broad profile-error reduction and failed target confirmation. This second experiment replaces the population generator rather than applying another response residual.

## Hypothesis and difference from prior experiments

Represent a perturbation as a change in the occupancy of real baseline cell states. Choose nonnegative weights on all query-fit control cells so their weighted aggregate approaches V63's predicted response. Then sample distinct original cells using those weights. All non-own-gene counts remain intact: no added Poisson noise, no moment-fitting edits and no global variance coefficient. The own-target gene alone is independently binomial-thinned at probability 0.1 before fitting and sampling. Thus libraries can decrease by the removed own transcripts; they are not refitted to original depths.

V59 estimated abundances in four learned states using donor treatments; V29 edited individual control counts. V68 instead optimizes individual-cell inclusion probabilities against the full nonpanel measured-gene response profile and leaves non-own counts untouched. No source treatment distributions or learned state labels enter this decoder. V63 remains the source-response predictor and is unchanged.

This is an independent adaptation of [stable balancing weights](https://doi.org/10.1080/01621459.2015.1023805) and [pivotal fixed-size sampling](https://doi.org/10.1093/biomet/85.1.89). Those statistical methods do not establish biological perturbation causality or competition performance. Our key assumption is restrictive: useful treated expression states must already be represented among controls. The method cannot create a response outside their convex hull.

## Optimization and sampling

Target probabilities are V63's closed desired bulk profile, with own suppression applied once. For each candidate cell, form a linearized TS50k bulk-profile discrepancy using its counts minus its library times the target probability. Minimize squared aggregate discrepancy on measured nonpanel genes plus a fixed quadratic penalty toward uniform weights. This is a convex surrogate, not direct optimization of the official expression score. The regularization coefficient is fixed at 0.0001.

Eight control-depth strata receive fixed sample quotas proportional to their original sizes. Each inclusion weight is nonnegative and at most 1/draw; each stratum has its assigned mass. The cap prevents collapsing onto fewer cells than the requested population. Accelerated projected gradient uses backtracking and up to 250 iterations; the convex duality gap and convergence flag are saved. Pivotal sampling preserves inclusion probabilities, selects the exact quota in each stratum and never duplicates a cell. The original depth groups are preserved; within-group depths and own-thinned final libraries can still vary.

Eight numerical tests verify the projection/cap/quota constraints, fixed sample sizes and inclusion frequencies over 6,000 draws, exact non-own count preservation, operator adjoints and finite-difference gradients, recovery of a known population shift, a zero-response uniform solution, and matched shuffled-weight distributions.

## Fixed biological test

Use all fit-control cells from the existing left partitions: H1, HCT116, HEK293T, CD4 resting donor 1 and CD4 stimulated donor 2. Check their profiles match the saved V63 query profiles exactly within numerical tolerance. Reference controls and all scored treatments are disjoint. The baseline response arrays retain their existing whole-study exclusion. These public datasets have been used in earlier experiments; no new independent study or blind validation is claimed.

Compare four arms on the first two unique target cases per context: unchanged V63, optimized real-cell sampling, equally concentrated weights shuffled within depth groups, and uniform real-cell sampling. The last two receive exactly the same own-target edit. Generate the same cell count as the old comparator (200 or 400); score identical treatment populations and independent real references. V63's ten saved comparisons may be reused only after checking hashes and bitwise truth equality. There are 40 treatment comparison records, including those ten reused baselines.

Also compare unchanged V63 and uniform real-cell sampling on a 400-cell zero-response null in each context: ten null comparisons. Actual counts are scored even when the convex-hull fit cannot achieve the desired mean, because the hypothesis concerns distribution fidelity. All raw Wilcoxon p-values, q-values, fold changes and profiles are retained for new count outcomes.

The primary must reduce paired equal-study profile error by at least 1%, satisfy fixed PDS/fidelity/reach/overlap/NMAE guards, improve at least two DE metrics, and beat both controls on profile error and at least two available DE metrics. Natural null calls must be no greater than max(5,V63 null calls) in each context. Undefined metrics stay missing and cannot satisfy a guard. The immutable plan records all thresholds. Two-target panels are coarse discrimination diagnostics, not reliable official-score forecasts.

Only passing development permits remaining canonical targets, a second technical seed, and then the independent reserve before packaging. A failed candidate is documented; no control arm, new regularizer or count mixture is promoted using these outcomes. D3 stays unopened in this initial test. No automatic upload or Git push.

## Reproduce

```cmd
cd /d L:\promec\Animesh\Download\Virtual
set "SCIENCE_PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
"%SCIENCE_PY%" -X utf8 test_v68_balancing.py
"%SCIENCE_PY%" -X utf8 v68_population_balancing.py freeze
"%SCIENCE_PY%" -X utf8 run_v68_experiment.py prepare
"%SCIENCE_PY%" -X utf8 run_v68_experiment.py counts
"%SCIENCE_PY%" -X utf8 run_v68_experiment.py nulls
"%SCIENCE_PY%" -X utf8 run_v68_experiment.py assess
```

Frozen scripts and compact receipts are versioned in the root and experiment analysis directory. Large matrices/count arrays are regenerable and ignored. Completed stages verify their stored hashes on resume. A separate result report records what completed and why any release was accepted or rejected.
