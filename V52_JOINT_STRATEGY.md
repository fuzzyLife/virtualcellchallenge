# V52: joint mean/rank response predictor

This experiment trains a shared neural response model from existing raw public perturbation populations. It follows the V51 diagnosis; it does not create or upload a competition submission. V24 remains the published reference until a new complete count-based evaluation establishes a useful candidate. D3 is not accessed.

The frozen executable specification is `experiments/v52_joint_mean_rank/analysis/plan.json`. All scientific code and source provenance used in the plan are hashed before label extraction or fitting. Later amendments, if needed for an execution failure, must preserve the original source and plan and explain whether any outcomes were known.

## Labels and data

Sources: H1, HCT116, HEK293T, CD4 Rest donor D2 and CD4 Stim48 donor D2. Existing extraction loaders verify original cached counts and cell metadata. HCT116 and HEK293T are one Orion study; both CD4 conditions are one CD4 study. These public studies have been repeatedly examined in earlier project work; the newly computed rank labels are not new blind biology.

For every supported competition target, sample two nonoverlapping treated populations and two nonoverlapping control populations within each source lane. Use at most 64 treated cells and 128 controls per lane per view, require at least 16 treated cells in each pooled view and at least 16 controls per lane per view. The sampled cell identities and lane weights are saved. Treated and control views are disjoint within each target; reusing controls across targets is allowed and is not counted as independent biological replication.

Compute full-axis labels on all measured output genes, excluding the intervention's own gene. Mean labels are batch-matched CPM log2 ratios, clipped to +/-3, using 0.1 only for the mean ratio. Rank labels are theta = U/(n_treated*n_control) - 0.5, with half credit for ties and no rank pseudocount. Average lane-specific statistics with treated-cell proportions. Both control-view means must exceed 5 CPM for a label to be eligible. Compute theta in every measured lane coordinate before this pooled filter, including coordinates that are rare in one lane. Two-view disagreements supply bounded reliability weights.

The source query profile is pooled source controls; mean labels use target-matched lane controls. These are related but not identical baselines. This residual source batch/domain mismatch is a limitation, not recipient treatment information available at inference.

## Model and exclusions

A shared 128/64 SiLU network with LayerNorm and 0.05 dropout predicts two bounded residuals: a mean log2FC residual around excluded-study V24, and a rank residual around an available-family rank pool. Mean residuals are bounded to +/-2 and final mean to +/-3; rank residuals to +/-0.25 and final theta to +/-0.5. Targets are not inferred from the own-gene suppression value.

Inputs are source-control basal expression/ranks, existing experimental donor mean effects and uncertainty, fixed STRING gene/target features, V24 mean predictions, and genuine donor rank effects/disagreement. All donor effects exclude the whole outer evaluation study and the whole fitting recipient study. Each raw donor study receives one rank vote; the two Orion contexts or CD4 conditions are not treated as independent study votes. Missing values have explicit indicators. The classifier/activity labels from V38 are not used as rank labels.

Use a target-name hash to exclude all occurrences of 20% of targets from fitting and retain them for source validation. Other allowed studies' experimental responses to those targets remain retrieval inputs, making this cross-context evaluation with potentially observed interventions, not unseen-intervention extrapolation.

Compare the joint network with a mean-only network, rank-only network, two independent linear residual heads, a trained target-label-shuffled joint network, V24 mean plus pooled rank, and a zero-rank baseline. The shuffled model uses a within-source cyclic permutation restricted to fitting targets, preserving gene alignment; unavailable shuffled labels abstain. Validation always uses genuine labels.

Train eight epochs, 192 batches per epoch, 2,048 pairs per batch, AdamW with learning rate 0.0005 and weight decay 0.001. All models use the same fitting-pair sample (at most 768 uniformly sampled eligible output genes per source target) and train-only feature normalization. Validation retains all hash0 fitting-source pairs and chooses only the checkpoint epoch. Whole-study test outcomes cannot select the epoch.

## Evaluation and scope

Three complete outer exclusions: H1, Orion, CD4. Score every eligible measured output gene for every supported target in the held-out sources, with explicit unique intervention counts. Report raw mean/rank MSE and MAE, strong-effect directional diagnostics, view noise and equal-study aggregation. These are not official FID/REACH or private anchor-scaled scores. Evaluate linear and shuffled controls under the same holdouts.

The prespecified pilot checks require joint mean MSE no worse than V24, rank MSE at least 5% below pooled rank, rank MSE below zero rank, and both errors below the trained shuffle. Passing is insufficient for packaging: significance, ties, row normalization, integer counts and official-style population metrics would still need separate decoder tests. A failed joint model remains a documented negative result, not a new submission.

After the outer tests, fit an all-source joint model using the same configuration and source-validation epoch rule; this is an exploratory deployment artifact even if the pilot fails. Its all-source validation is not an independent test. Models, normalization and full-axis held-study mean/rank arrays are saved with hashes; large binaries stay outside Git.

## Reproduce

Use the pinned Python 3.12 scientific runtime and requirements-v35-science.txt. SCIENCE_PYTHON is its interpreter path. Keep the official Python 3.14 CLI environment separate.

~~~bat
SCIENCE_PYTHON test_v52_joint_response.py
SCIENCE_PYTHON run_v52_joint_response.py freeze
SCIENCE_PYTHON run_v52_joint_response.py prepare --workers 2
SCIENCE_PYTHON run_v52_joint_response.py train --outer H1
SCIENCE_PYTHON run_v52_joint_response.py train --outer Orion
SCIENCE_PYTHON run_v52_joint_response.py train --outer CD4
SCIENCE_PYTHON run_v52_joint_response.py evaluate
SCIENCE_PYTHON run_v52_joint_response.py train --outer all
~~~

Freeze only once in an untouched experiment directory. Existing labels and model files resume only when their checksums match. Do not overwrite an existing frozen plan or silently revise a failed experiment.
