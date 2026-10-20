# V76 results: shrinking the complete response did not generalize

V76 tested normalized mixtures of the full V66 response, including own-target knockdown, with the unperturbed profile. The development choice was mixture075, selected with both expression and discrimination checks. The original expression-only ranking is retained in `mse_only_preselection.json`; no candidate was changed after confirmation.

On the broader 655-pair check, desired-profile error ratio was 0.982708391 but discrimination fell by 0.026819978. On all30 fixed count cases, PDS fell from0.759126984 to0.688690476. On the15 remaining-target profile/DE cases, fidelity improved from0.568444162 to0.575940266, but reach fell from0.111707816 to0.069288284 and Jaccard from0.084688087 to0.081329376. The paired profile-error ratio was0.945719588. Full three-study NMAE was unavailable; missing values were not treated as zero.

All90 baseline/candidate/shuffled count records passed independent replay of four DE statistics, BH arrays, PDS, and720 raw-rank coordinates. Every regenerated V66 baseline matched the original frozen cache. These checks establish numerical integrity, not improved predictive performance.

The critical guards failed. No repeat seed, production fit, upload file or replacement candidate was released. H1/HCT116 have only four unique cases each; V75's three-development/one-confirmation split could not support held-out PDS. V76 used all30 cases for PDS and the remaining15 for primary expression/DE, with the reuse explicitly stated. V78 onwards use two development targets per context, leaving at least two for confirmation.

Reproduce saved statistics with the scientific Python environment:

```cmd
python -X utf8 audit_v76_counts.py
```

The model, fixed protocol and selection are in `v76_response_dilution.py`, `V76_STRATEGY.md` and `experiments/v76_response_dilution/analysis`. Build/package adapters are gated preparation code and were not executed.
