# V72 results: coverage fixed, predictive transfer rejected

The corrected batch-offset estimator restored supported source targets from 12 to 25 in H1, 13 to 275 in HCT116, and 23 to 287 in HEK293T. Treated-cell eligibility and counts were unchanged. Synthetic small-batch tests and the complete source-exclusion/correction audit passed. Coverage details and hashes: `coverage_comparison.json`.

The fixed development screen nevertheless rejected every candidate. Profile-error ratios versus V66: raw25 **1.032179512** (+3.22%), center25 **1.030686921** (+3.07%), center50 **1.100953584** (+10.10%), credible25 **0.999906708** (-0.0093%). The credible change is below the frozen 0.5% gain criterion. Raw25 and center25 also lost to their matched label shuffles. No candidate was selected, no confirmation case opened, and no V72 production artifact generated.

This distinguishes data coverage from useful transfer: more genuine source evidence alone did not fix context/representation mismatch. V73 adds a prior centered on unit response for weak/high-dispersion count likelihood and transfers effects in bulk assay space, with a raw-scale comparator. It is a new frozen experiment, not a revision or promotion of a failed V72 arm. All biology is reused public data, not blind independent validation. V66 remains best published.
