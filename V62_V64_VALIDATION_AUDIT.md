# V62-V64 validation audit: case counts, target identity and repeated biology

The user requested deeper investigation and a new submission to test. Three new response representations, a target-margin guard, and scalar-only calibration were implemented. None passed the preset public count acceptance rules. Their gains and failures are preserved; a new test artifact must therefore be explicitly experimental, with V24 retained as the best published reference.

## A concrete validation gap

The shared V34/V37 cohort has **44 predefined case configurations**, but **30 unique target/context combinations** and **34 distinct actual-cell row sets**. Different row sets may overlap and are not biological replicates. Repeating a generator seed does not add biological evidence.

| Context | Configurations | Unique targets | Distinct measured row sets | Actual cells/configuration |
| --- | ---: | ---: | ---: | ---: |
| H1 | 8 | 4 | 8 | 400 |
| HCT116 | 8 | 4 | 4 | 400 |
| HEK293T | 12 | 6 | 6 | 400 |
| CD4 Rest D1 | 8 | 8 | 8 | 200 |
| CD4 Stim48 D2 | 8 | 8 | 8 | 400 |

Some HEK293T additional panels contain the same perturbation label twice. Treating those as separate competitors does not match the competition's unique-target PDS definition. The new auditor keeps the first case_index per target within each context/seed panel and saves corrected panels separately. Original count arrays and original failed verdicts are unchanged. Development panels already have four unique targets, so this correction cannot erase their failures.

For uncertainty and mean-error aggregation, first average repeated case/seed errors within each target, then compare summed errors within a context, then average contexts within each of the three studies. The revised count-profile summary records all replication levels. Earlier pre-deduplication script/summary/audit values are retained in pre_deduplication_analysis.json. Nothing here creates an independent new validation set.

## What the experiments say

V62 source-common removal alone gives a smaller mean benefit than the assay-space estimator. The assay estimator transforms uncertainties but reuses an LFC prior grid and variance floor, so its gain cannot be attributed to a units-invariant coordinate transformation. V63 adds a fixed prior-template margin safeguard. V64 preserves V24 gene signs and locations and changes only vector magnitude; its mean-side promise fails the final-cell test. This rejects a blanket claim that the more elaborate method or a simple scalar automatically improves competition performance.

The actual V63 integer-count comparisons reduce paired plug-in expression-profile error across every public context, but lose discrimination in the small count panels and do not pass all overlap/reach checks. Its broader desired-mean identity loss is smaller than V62's. These are reasons to study a controlled test artifact, not evidence of an improved hidden score.

Weak perturbation-to-control distances make raw expression-error ratios unstable, especially HCT116. The added target-weighted error analysis instead divides by V24 prediction error for a paired diagnostic. It is not the official jackknife/anchor score. Means, realized cells, target identity, significance calls and count variance remain separate diagnostics. No official overall score or ranking is forecast.

## What remains missing

We still lack a genuinely independent, adequately powered new-study benchmark for these revisions. All five public contexts and their labels have been reused; HCT116/HEK293T constitute one study and the CD4 conditions one study. The real challenge asks for unseen-context intervention responses; better within-study reconstruction or control matching does not identify those responses. D3 remains unopened. The inherited count generator's excess null DE is also unresolved.

Relevant root scripts are audit_v64_unique_panels.py, audit_v64_unique_mean_control.py and analyze_v62_v63_diagnostics.py. Corrected results preserve canonical target IDs and all measured zeros. The research and model protocols are V62_CONSERVED_RESPONSE_STRATEGY.md and V63_IDENTITY_GUARD_STRATEGY.md; V64_RESULTS.md documents the rejected scalar follow-up. The main reproducibility receipts are under experiments/v63_identity_guard/analysis and experiments/v64_response_magnitude/analysis.
