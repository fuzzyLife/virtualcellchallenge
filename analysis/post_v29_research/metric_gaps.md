# Leaderboard metric gap from the supplied screenshot

Reviewed 2026-10-15. Scaled components are higher-is-better. Competitor values are manually transcribed rounded screenshot values; their anchor IDs are not visible. These are descriptive comparisons.

| Component | V24 | iCell | X-89 | iCell minus V24 |
|---|---:|---:|---:|---:|
| Perturbation discrimination | 0.6442 | 0.874 | 0.776 | +0.2298 |
| Expression accuracy | 0.0476 | 0.787 | 0.458 | +0.7394 |
| DE fold-change accuracy | 0.0823 | 0.419 | 0.283 | +0.3367 |
| Direction fidelity | -0.0158 | 0.153 | 0.045 | +0.1688 |
| Direction reach | 0.0848 | 0.444 | 0.362 | +0.3592 |
| Significance overlap | 0.0149 | 0.088 | 0.088 | +0.0731 |

V24 overall: 0.142997855. Screenshot iCell: 0.4610; X-89: 0.3353.

Expression and fold-change accuracy account arithmetically for 56.4% of the rounded six-component gap to iCell. This does not identify its method or predict the gain from changing one part of our pipeline.

The exact V24 raw expression error is 0.960415; its raw DE fold-change error is 0.950686. These raw errors are lower-is-better.
