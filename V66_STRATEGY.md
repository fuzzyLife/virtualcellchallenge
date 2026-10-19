# V66: strength-matched original assay response control

V65's response screen rejected both half-strength normalized variants but selected its three-quarter variant, with only 0.27% lower desired-profile error than V63. Strength and normalization therefore remain confounded. V66 tests the original successful assay representation at the same correction strength, retaining V24's 80% prior-margin guard and original residual10 count generator. It does not add a new ensemble or claim normalized data are more accurate.

The fixed candidate expands the unguarded V62 assay_center50 correction from V24 by 1.5 before applying the unchanged guard. V63 at the original strength must replay bitwise in each context. A cyclic source-target derangement receives the same expansion and guard. Whole H1, Orion (both HCT116 and HEK293T), and CD4 studies are excluded from their respective response banks. Query controls are legitimate inputs; query treatments cannot enter predictions. The adversarial boundary audit verifies source inputs remain identical when recipient outcomes/support/uncertainty/knockdown/profiles are poisoned.

This outcome-informed follow-up was frozen after V65 means and after V65 count jobs began. At specification time, only runtime/progress logs had been inspected; no new count accuracy summary was used. These are repeatedly reused public datasets, not blind biological validation. The V65 unique-target test machinery is reused without numerical changes: 30 unique target/context pairs, five contexts and two technical generation seeds, four arms including V24 and V63. D3 remains unopened.

The cross-version choice is frozen in V66 plan.json before V66 count outcomes. Prefer a candidate passing every original count guard; if both pass, choose the smaller equal-study generated-profile error ratio. If neither passes, an explicitly experimental diagnostic is eligible only with no study profile error >3% worse and at least two improved raw metrics versus V63. No post-count promotion of another strength. `select_v65_v66.py` applies the rule; all failures must be reported. Neither local proxy nor geometry preservation guarantees a private leaderboard gain.

## Reproduce

Use the same documented scientific Python 3.12 environment as V65; `.competition_env` is for official packaging.

```cmd
SCIENCE_PYTHON v66_assay_strength.py freeze
SCIENCE_PYTHON v66_assay_strength.py prepare
SCIENCE_PYTHON audit_v66_boundaries.py
SCIENCE_PYTHON run_v66_experiment.py means
REM If the fixed mean screen selects assay75:
SCIENCE_PYTHON fit_v66_submission.py freeze
SCIENCE_PYTHON run_v66_experiment.py parallel
SCIENCE_PYTHON audit_v65_results.py
SCIENCE_PYTHON audit_v66_results.py
SCIENCE_PYTHON select_v65_v66.py
REM Fit/generate/package only the selected version, explicitly experimental if needed.
```

`finish_v66_screen.py` sequences these computations and waits for V65 counts to complete before launching its two count workers. The separate production adapters preserve full response replay, packaged count/axis checks, original-control identities/guides/libraries and 12 full-group regeneration checks. Nothing uploads or pushes Git automatically. See the complete [18-submission review](V65_SUBMISSION_HISTORY.md), [V65 protocol](V65_STRATEGY.md), frozen plans and compact receipts for exact source/file hashes.
