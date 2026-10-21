# V86: covariance-mixture posterior transfer

The new hypothesis is to learn how response effects co-vary across source contexts, then condition on the measured effects available for each target/output coordinate. A Gaussian-process context kernel supplies uncertainty for an anonymous query context using only NTC expression. This is different from averaging predictions, fitting another response amplitude, or V84's learned mean operator on a nonoverlapping essential-target panel.

This is an independent adaptation inspired by [multivariate adaptive shrinkage](https://www.nature.com/articles/s41588-018-0268-8) and its [author documentation](https://stephenslab.github.io/mashr/). We did not execute mashr or reproduce its association-study benchmarks; those results do not establish CRISPR or competition accuracy.

## Inputs and exclusions

Use the existing measured K562, CD4, H1, HCT116 and HEK293T response cache, with its measurement masks and uncertainties. Apply the same bounded source dose normalization used by V66. Convert responses into the competition's log1p(50000 times bulk proportion) assay units, propagating uncertainty through the derivative and retaining a 0.0025 variance floor.

For each recipient study, exclude all its response data before fitting **and** conditioning: H1; Orion comprising both HCT116 and HEK293T; or CD4 comprising every state and donor. All 300 panel genes are excluded as model outputs. Missing coordinates are absent observations, not zero effects.

The source-context feature axis contains nonpanel genes with finite control expression above 5 CPM in every allowed source. Center and normalize their rank vectors. Query features rank only finite measured values and put zero in missing **centered rank embedding coordinates**, not in an RNA matrix. At least 90% of the frozen source feature axis must be measured or the model abstains. No pseudocount is added for ranking.

## Model and fixed controls

A source covariance component is `scale^2 * (0.8 * kernel + 0.2 * I)`. Kernels are Gaussian functions of rank-feature distances with length 0.05, 0.2, 1 or 4, plus an all-ones shared kernel. Scales are 0.025, 0.1, 0.4 and 1 in assay units. A null component gives 21 components in total. The query retains its own unobserved 0.2 idiosyncratic variance; apparent source similarity does not make the query response known.

Fit mixture weights to 20,000 deterministic source target/output coordinates with at least two valid sources. Each Gaussian marginal likelihood uses its actual observed-source submatrix and sampling variances. Run at most 100 mixture EM updates, with a fixed 0.5 pseudocount per component and a 1e-7 maximum-weight-change stopping criterion. Report the iteration limit if reached. Correlated genes and targets are working-likelihood observations, not independent biological replicates; posterior intervals are modeling quantities without established coverage in the target contexts. The code uses mean +/-1.2815515655 standard deviations, not exact mixture quantiles; fitted mixture weights are plug-in estimates.

The primary `interval25` moves the V66 assay response toward a nominal 80% interval computed from the conditional mixture mean and variance using a normal approximation, converts it back to fold change, and applies one quarter of a correction bounded at +/-0.5 log2 units. `posterior25` instead uses the posterior mean. `shared25` independently fits weights restricted to the shared/null covariance components, providing a strong simpler comparator. `interval25_shuffled` refits the entire prior after independently shifting source target labels and conditions on those shifted responses. This destroys target correspondence while preserving each source's marginal distribution; missing-support patterns can change.

All use the existing 80% V66 identity-margin guard and the same count decoder. Own/panel effects remain unchanged. Unsupported source coordinates abstain. No arm or strength is selected after recipient results; `interval25` remains the fixed primary.

## Tests and decision

Synthetic tests verify Gaussian conditioning against direct linear solves, missing-source poisoning invariance, exact zero-response means, and transfer to an unobserved synthetic context. The latter improves MSE versus copying a source; this tests the implementation, not real-world performance.

Development includes all 655 reused public pairs. Require at least 0.5% lower equal-study profile error, no study more than 3% worse, PDS loss at most 0.005, and lower profile error than both shared-only and refitted-shuffled controls. Whole-study exclusion is retained; these are reused public biological data.

The count diagnostic always tests the fixed primary on the existing 14-target D3 cohort with matched 400-cell actual/template populations, 5,006 inference NTCs and a separate 4,994-cell reference. D3 was used in V84 and is not blind. The entire CD4 response family remains excluded from V86.

The original 14 V84 V66 baseline populations are reused only after exact equality checks on V66 effects, source CPM, bulk proportions, gene/target axes, cohort and seed, and unchanged emitter/metric code. Their current truth/reference profiles are replayed. This avoids claiming that identical regenerated baseline cells would be new evidence. Emit 42 new populations: primary, shared-only and refitted-shuffled predictions for each target. Independently replay profiles, BH, DE and panel discrimination for all 56 records. There are 336 new raw tied-rank coordinate checks plus 112 previously checked baseline coordinates.

Require profile error not to increase, PDS loss <=0.005, fidelity loss <=0.002, reach loss <=0.005, Jaccard loss <=0.0005, NMAE <=1.01 times baseline, at least three improved endpoints, and better profile error than both controls. Run seed offset 620000 only if both development and the first count test pass. No upload package otherwise; no raw public metric is converted to a fabricated hosted overall score.

## Reproduce

Use the existing Python 3.12 scientific runtime, not the official packaging environment. Plans, code hashes and source-cache hash are fixed before fitting. Existing fits and prediction files refuse overwrite.

Windows CMD (matching the submission terminal):

```bat
cd /d L:\promec\Animesh\Download\Virtual
set "SCIENCE_PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
"%SCIENCE_PY%" -X utf8 test_v86_posterior.py
"%SCIENCE_PY%" -X utf8 finish_v86_experiment.py
```

`requirements_v85_v86.txt` records the tested CPU package versions. A replacement Python3.12 environment can install these requirements; point `SCIENCE_PY` at that environment. The official VCC packaging environment is separate. Required prior data/artifacts must be present and match their recorded hashes.

The runner fits genuine and shuffled models for each of three outer studies with two workers, predicts all public contexts and D3, reports all development arms, audits the model and boundaries, and runs the paired count test. Individual stages are in `run_v86_experiment.py` and `run_v86_counts.py`. `audit_v86_experiment.py models` checks source likelihoods against SciPy, conditional posteriors by direct solves, positive semidefiniteness, recipient-study poisoning and realized identity margins. `audit_v86_counts.py` replays the final count statistics with arrays cached in memory.

The committed final plan is sufficient to start in a fresh analysis checkout with the required prior data and artifacts. `freeze` is a once-only stage used when creating a genuinely new experiment directory. Large arrays, model weights and logs are ignored by Git; compact model metadata, provenance, failures and audit receipts are retained.
