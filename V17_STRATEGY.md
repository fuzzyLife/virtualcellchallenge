# V17: learning intervention-specific population variation

V17 tests whether a learned distribution of response strength improves on V16 while keeping its average perturbation predictions fixed. It is an exploratory transfer experiment on public source biology. A more spread-out scatter plot alone is not evidence of better predictions.

## Literature and the resulting hypothesis

- [PopPert, September 2026 preprint](https://arxiv.org/abs/2609.01357) models population distributions with gene marginals and shared low-rank dependence. Its useful lesson here is to evaluate joint distributions rather than only average profiles. We did not copy its code or implement its full neural mixture/copula model.
- [Song et al., Nature Cell Biology 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11906366/) estimates a perturbation-response score per cell, motivating a shared latent response strength across affected genes. Its genome-scale Jurkat experiment has multiple guides per cell; those measurements require joint-guide modelling before reuse as single-target responses. Our cross-fitted variance estimator is an independent, simpler adaptation, not the published PS implementation.
- [CausalPerturb author code](https://github.com/zcfu21/CausalPerturb) separates basal state and heterogeneous treatment effects. Its fitted cell effects require observed treated cells; they cannot identify hidden A/B/C treatment responses from controls alone.
- [CellOT](https://www.nature.com/articles/s41592-023-01969-x) motivates distributional comparisons of unpaired populations. V17 is not a trained CellOT transport network.

We hypothesize that some targets produce a mixture of weaker and stronger responders. Sampling a shared response gain per cell creates coordinated changes across response genes. Independent per-gene noise cannot represent that structure.

## Public data checked

Root `collect_v17_resources.py` pins Arc PIE essential-screen metadata; `collect_v17_public_tables.py` collects author manifests, metadata and target coverage. URLs, release revisions, MD5/SHA-256 checks and file sizes are saved in `experiments/v17_learned_heterogeneity/analysis/`.

| Resource | Measured overlap with our 300 targets | Current use |
|---|---:|---|
| Arc PIE Replogle/Nadig essential screens, metadata revision `91f86154aa1761a2dd0d885b979b9169f3b9983f` | 0 of 1,857 perturbation names | Metadata audit only; no new matched response examples |
| [Multi-donor iPSC CRISPRi count release](https://doi.org/10.6084/m9.figshare.27989294), fitness / non-fitness screens | 31 / 166 | Downloaded cell metadata; medians 50 / 22 assigned cells per target pooled across donors; no target reaches 400 cells |
| Same iPSC targeted screen | 5 | EZH2, GIGYF2, NT5DC1, STT3A, ZNF526; 1,133–1,923 cells pooled, only 96–144 maximum in one donor |
| [KOLF2.1J atlas](https://www.nature.com/articles/s41587-026-03199-w), author browser | 273 targets, 18,106 output genes | Pinned summary metadata and author export code; summaries are quantized mean NTC z-scores, not counts or log fold changes |

The KOLF full QC-filtered AnnData is 189,393,177,972 bytes. We did not download it. The multi-donor iPSC RNA matrices are also not downloaded or fitted in V17. Unassigned guides are not controls. Donors, biological batches, assay differences and response magnitudes must be preserved before either resource becomes a production source. These newly collected metadata are separated from the frozen V17 tests; they are not claimed as new cell-level training examples.

V17 learns from the existing verified public raw-count caches: H1, HCT116 and HEK293T. K562/CD4 currently supply mean-response evidence, not cell-level variation training.

## Model and fitting steps

1. Verify original sparse count matrices, selected cell identities and gene-axis masks using `fit_v7_response.raw_inputs`. Normalize on the complete available challenge-aligned library and apply log1p CPM for the variation estimator. Missing output genes remain masked.
2. Within each source, randomly split controls from every batch into two halves. Batches with fewer than eight controls in either half are excluded. Require at least 80 treated cells per target after this eligibility filter.
3. Split each target's treated cells into two halves. Estimate its response direction from one treated half and one matched-control half. Subtract the control mean of each cell's own batch before pooling. Exclude all 300 competition panel genes; select at most 256 standardized response genes above the fixed signal threshold, with at least ten genes.
4. Project the other treated half and its independent control half onto that direction. Estimate excess treatment variance above matched within-batch control variance, correcting for uncertainty from estimating the control mean. Divide by squared mean response, cap CV squared at 0.5, shrink for sample size and signal-to-noise, and abstain below SNR four. Swap halves and average.
5. Average target estimates equally over available training families. Every source-family benchmark excludes the entire test family from both response and variation fitting. Missing or unsupported variation is exactly zero, which returns V16.
6. At generation, draw one bounded positive mean-one gain per cell. Multiply the target's response vector by this gain. Keep V16's Poisson coefficient 0.65, empirical coefficient 0.10, own-target factor 0.10 and sampled library sizes fixed. Dual-moment fitting retains the original V16 target population mean and bulk constraints. Integer counts and libraries are checked.

The production fit has nonzero variation for 85 targets. This models a shared latent response strength; it does not infer causal dependence on a cell's basal state or learn every possible covariance change. Survival bias and assay differences can affect source variation estimates.

## Prospective tests and release policy

`prepare_v17_experiment.py` freezes the development identities, candidates, dependencies and gates. Fourteen target/family pairs, two seeds and 400 real/predicted cells produce 140 comparisons: V16, half/full learned variation, generic variance and target-shuffled variance. Control fitting and evaluation references are disjoint. Six real NTC cases check bitwise V16 identity. Generator invariants also check zero-variance identity, determinism, integer counts, moment constraints and exact libraries.

Keep raw PDS, fidelity, reach, Jaccard, NMAE and expression error separate. Apply existing available-metric and worst-family guards, require three improvements including DE, require improvement in response-projection Wasserstein distance, and require superiority to generic and shuffled controls on at least two DE metrics. Missing NMAE is disclosed, never filled with zero. Raw local metrics are not official anchor-scaled leaderboard scores.

Only a passing finalist opens reserved algorithm verification. `prepare_v17_confirmation.py` freezes eight target/family pairs disjoint from current development, selected alphabetically among groups with at least 400 cells. The original draft accidentally reused three development pairs; the corrected verification plan was saved before candidate selection. HCT116 has only two eligible unused names, CAT/MYRF, already inspected in earlier versions. This is public algorithm verification, not blind biological validation. The original development plan and all development inputs remain unchanged.

The first input-power pilot required too many cells within individual cached batches and produced no eligible H1/HCT116 estimates. Before candidate benchmarking it was replaced by batch-centred pooling with independent control halves. A subsequent eligibility fix skips underpowered HEK293T batches instead of aborting the source. Original plans/code and amendment explanations are preserved in `analysis/input_power_pilot/`.

## Reproduce from Windows CMD

Run from `L:\promec\Animesh\Download\Virtual`. Commands creating plans/models refuse to overwrite them; use a fresh experiment directory for a genuinely new experiment.

```bat
.\.competition_env\Scripts\python.exe collect_v17_resources.py
.\.competition_env\Scripts\python.exe collect_v17_public_tables.py
.\.competition_env\Scripts\python.exe test_v17_heterogeneity.py
.\.competition_env\Scripts\python.exe prepare_v17_experiment.py
.\.competition_env\Scripts\python.exe run_v17_experiment.py models
.\.competition_env\Scripts\python.exe run_v17_experiment.py null
.\.competition_env\Scripts\python.exe run_v17_experiment.py development --workers 6
.\.competition_env\Scripts\python.exe prepare_v17_confirmation.py
.\.competition_env\Scripts\python.exe release_v17_submission.py
.\.competition_env\Scripts\python.exe plot_v17_evidence.py
```

The saved plan is the exact executed input specification. The correction history is retained; do not silently regenerate a different plan over it. The release runner evaluates the frozen finalist once, independently audits source exclusions and identities, and only then fits/builds/packages/validates production output if all gates pass. Production retains V16 control sampling seeds, permitting exact reuse of groups whose response and learned variation are unchanged. Full-payload comparison and deterministic replay must pass before an upload is recommended. No automatic upload is performed.
