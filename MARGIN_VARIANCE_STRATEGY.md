# V78-V80: preserve the published response while testing cell distributions

V69's published regression showed that frozen intended effects do not protect realized count profiles. These experiments instead take the exact regenerated V66 count population as the reference, including its original row libraries and gene sums. This is a controlled decoder experiment; there is no new mean-response learner.

## Invariants and limitation

The V29 projection first fits probabilities, integerizes with exact row and column totals, then uses integer checkerboard swaps to repair the cell-normalized mean. Every returned population must have exactly the same row libraries and gene totals as its V66 parent. Raw bulk expression profiles and PDS are therefore exactly protected for any true response, not merely approximately fitted on public cases. The cell-normalized mean must have L1 difference at most0.0005, otherwise the parent is returned. Movement limits and fallback counts are recorded.

Official expression accuracy includes a sampling correction; changing cell distributions can change that correction even at identical bulk profiles. Differential-expression means, raw rank p values, FID, JAC and REACH may also change. These quantities require tests. Exact bulk identity is not a guarantee of identical official MSE or improved leaderboard score.

## Three distinct proposals

- V78: control-based Gamma-Poisson shape coupling at10%/25%, a study-excluded learned shape at25%, and stronger50%/100% query-shape proposals. All300 panel genes are excluded from the direct proposal. Smaller arms cap alpha changes; full-query arms do not. Existing V69 shape machinery is read only; no sealed V69 artifact is modified.
- V79: integer molecule thinning plus negative-binomial innovations. For retention rho and signed parent excess variance a0, innovation excess variance is `(a_target-rho^2*a0)/(1-rho)^2`. Unlike a positive NB floor, this can account for an underdispersed parent. The algebra and stochastic moments are tested independently. Only the current own gene is excluded from the proposal.
- V80: fit control detection probability against library size and positive excess-count dispersion. Transport the empirical zero probability by a relative NB detection response, then generate positive counts. Compare control-cell-coupled detection, independent detection, and empirically scaled whole-control counts. This addresses the fact that matching two moments alone does not specify the ranks or zero fractions used by DE scoring.

Direct-proposal exclusion does not mean a gene's individual cell counts are fixed: joint repair may use carrier genes. All original row/column invariants remain enforced.

The code inspection motivating these tests found that V66's initializer is25% expectation,65% Poisson draw and10% empirical control counts. That blend attenuates biological variation; even its conditional sampling variance can be below Poisson. This is a plausible modeling concern, not proof that restoring full control variance will improve prediction.

## Prospective evaluation

All source treatment responses exclude the complete recipient study: H1, Orion (both HCT116 and HEK293T), or CD4 (all donors/conditions). Query shape/hurdle parameters use only the query-fit NTC controls. No query treatment or held-out reference controls are used to fit them. Query-fit NTC sample sizes are768,1742,3558,5007,5005 across the five contexts; challenge controls have about18400 per context. That difference limits transfer claims.

The first two unique targets per context are development (10 populations). The remaining2/2/4/6/6 are confirmation (20); every PDS panel therefore has at least two targets. Confirmation is conditional on a selection passing the original frozen development criteria. Two planned seeds are0 and620000; they are technical repeats, not independent biological replications. Public biology has been used repeatedly and is not blind.

V78 requires FID gain0.002, JAC gain0.0005, reach decrease at most0.005, NMAE on explicitly available studies no worse than1%, and genuine-shape superiority to shuffled controls. V79/V80 separately prescribe a precision-focused FID gain0.01, JAC gain0.001 and reach decrease at most0.02, with the same NMAE, exact-profile and shuffle checks. Each hypothesis and its different objective was frozen before its own outcomes; failed candidates are not promoted by changing these cutoffs.

Every baseline is regenerated and compared against the original V66 arrays. Each generated case retains raw p, adjusted p, effect sizes, eligibility, profiles, jackknife variances and truth. Eight eligible coordinates per population are checked against an independently computed tied-rank Mann-Whitney statistic. `audit_shape_trials.py` independently replays BH, all four DE metrics and PDS for the complete saved arrays and verifies exact-profile identity.

## Reproduction

Use the verified scientific environment (Python3.12.10, NumPy2.5.4, SciPy1.18.1, pandas3.0.6, h5py3.16.0, scikit-learn1.9.1, Numba), not the separate official CLI environment. Existing plans and selections are immutable; freeze commands deliberately refuse to overwrite them. Saved statistical replay is safe:

```cmd
python -X utf8 test_v78_margins.py
python -X utf8 test_v79_innovations.py
python -X utf8 test_v80_hurdle.py
python -X utf8 audit_shape_trials.py --version 78 --stage development
python -X utf8 audit_shape_trials.py --version 79 --stage development
python -X utf8 audit_shape_trials.py --version 80 --stage development
```

`run_v78_development.py`, `run_v79_development.py`, `run_v80_development.py` execute bounded-memory source jobs and the fixed selection. Their per-case counts action can resume completed arrays after checking hashes. Original query-control fits, dependency hashes, arm definitions and execution protocols are retained under each version's analysis folder. Large count arrays are local reproducible outputs excluded from Git; compact records and code are committed.

## V81: opposite-direction diagnostic

After inspecting completed V78 and partial V79/V80 outcomes, V81 prospectively tests the opposite direction, reducing count residuals instead of restoring query-control variance. The proposal is `(1-gamma)*X + gamma*library*bulk` at gamma0.25/0.50, with negative gamma (clipped at zero before projection) as an opposite-direction control. Exact V66 row/gene sums and the same CPM contract still apply. No new mean model or source coefficient is learned. This is explicitly decoder calibration, not a claim of biologically realistic low variability.

V81 has its own frozen development criteria: FID+0.002, JAC no worse0.0005, REACH no worse0.005, available NMAE no worse1%, and FID/JAC superior to opposite expansion. This is not promotion of a failed V78-V80 candidate. If selected, both20-case confirmations must pass, followed by a control-only false-DE guard. In each of five disjoint original NTC null splits, candidate false calls must be at most1.10 times the paired baseline plus10. This allowance was pinned before any V81 null outcome. Failed development skips confirmation/null/production; it is not evidence that null calibration passed.

Additional code: `v81_residual_temperature.py`, `run_v81_variance.py`, `run_v81_development.py`, `check_v81_nulls.py`. Frozen plans and source hashes are in `experiments/v81_residual_temperature/analysis`.
