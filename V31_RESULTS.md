# V31 guide-resolved response results

**Decision: no V31 release.** Guide-aware means failed the frozen screen. V24 remains the best published project entry at 0.142997855; no generated counts, fresh donor reserve or package were opened.

Whole-guide-held-out calibration within public D2 selected correction strength 0.25. There are 86 held-guide cases evaluated at three strengths (258 calibration rows), not 258 independent biological experiments. After minimum-cell filters, 22 resting and 21 stimulated targets support multiple eligible guides;34 targets are represented per condition.

## Mean-only development

|Arm|All-study LFC MSE ratio|D1 donor LFC MSE ratio|D1 donor CPM MSE ratio|
|---|---:|---:|---:|
|guide|0.996211297|0.988633890|0.964154925|
|pooled|0.996031583|0.988094749|0.962708612|
|guide_shuffled|0.995928393|0.987785179|0.960694679|
|target_shuffled|0.999949930|0.999849789|1.022395183|

Guide-aware transfer reduces D1 squared fold-change error 1.14% and mean-CPM squared error 3.58%. Ordinary pooled responses and shuffled guide labels do slightly better. It beats target-shuffled fingerprints, supporting some useful target-specific donor response information, but does not establish a benefit from separating guide effects. The all-study fold-change error improvement is only 0.38%.

H1 and Orion mean profiles remain exactly at baseline because the compatibility gate is 0. The D1 donor gate is 1. These are fallback/noninferiority checks on H1/Orion, not successful guide-response generalization to those lineages.

Failed checks: overall squared/absolute fold-change gains, the required 5% donor squared-error gain, and beating pooled/guide-shuffled controls. Worst-family noninferiority, mean-CPM noninferiority and beating target shuffle pass. Diagnostic pooled/shuffled arms are not promoted after the sole guide candidate fails.

The 36 development cases include 28 genuine 400-cell H1/Orion cases and eight 200-cell D1 donor cases. Four arms give 144 deterministic mean comparisons; **no synthetic cell populations were generated**. Repeated seeds reuse the same biological contrasts. These approximate clipped/smoothed mean errors are not official Wilcoxon/anchor-scaled scores.

## Validation boundary and verification

The raw teacher uses D2 resting/stimulated cells only. D1 validation excludes every CD4 summary from the baseline, then adds only raw-D2 guide-derived information. This is a donor-transfer test within an existing study, not a whole-CD4 or unseen-lineage holdout. D2 stimulated cells are not validation. The untouched D3 reserve was not accessed.

Three guide-pooling invariants pass. Independent dense raw-count calculations reconstruct selected source guide means, knockdown proxies, lane-control variances, hierarchical effects and confidence. Held-guide roles, selected strength and study-level arithmetic are recomputed; route/source-family exclusion and own-gene preservation are checked. Detailed measured tolerances are archived in the independent audit.

The wider conditional-transport and full transition fine-tuning proposals remain untrained. The staged plans require a useful mean finalist before those larger downstream runs; these rejected mean pilots do not demonstrate that those model classes fail.

Reproduction and exact assumptions: [V31_STRATEGY.md](V31_STRATEGY.md), root `v31_methods.py`, `run_v31_experiment.py`, `test_v31_guides.py`, `audit_v31_experiment.py`. Prior frozen models and V24 package are unchanged.
