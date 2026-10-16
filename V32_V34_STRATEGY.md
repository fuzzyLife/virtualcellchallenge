# V32–V34: three sequential model tests

These pilots test three distinct hypotheses in order: paired perturbation-response transfer (V32), a tabular foundation model (V33), and a conditional RNA flow that learns both response means and variation (V34). V24 remains the published comparator. No previous scientific files or packages are replaced, and no submission is uploaded.

## Shared validation boundary

H1, Orion (HCT116 and HEK293T together), and CD4 (every condition/donor) are held out as complete families. The public cases have already been examined in earlier experiments; these results are development evidence, not blind confirmation. K562/RPE1 broad pseudobulk sources are available to the paired/tabular fits. Their rows are fractional expression summaries, not raw single cells. The D3 donor reserve stays unopened.

Source-side validation uses a deterministic target-identity hash. Its 20% validation identities do not enter source response bases or fitting labels. A measured response for the query target in an allowed source is an admissible query token for V32/V33. No treated response from the recipient family is available at inference. The three approaches are all tested regardless of earlier failures.

## 1. V32: paired, context-conditioned response mapping

Learn directed changes between responses to the SAME targeted gene in different public contexts. Inputs include the observed source-response program scores, source/recipient control-rank coordinates, STRING target features, and explicit context interactions. Fit a shared 16-program basis on up to 2,048 training-selected, mutually measured genes. Use exact own-output-masked ridge. Compare paired labels with within-route shuffled identities, direct source copying, and V24.

A nested whole-family loop selects the V24 correction strength from 0.25, 0.5, and 1. Ridge regularization is selected using source-only identity-held examples. Require at least 1% squared-LFC and 0.5% absolute-LFC error gains, worst-family error <=1.02, and better paired performance than the shuffle and copying controls.

The initial execution repeatedly rehashed multi-gigabyte V24 source artifacts. `v32_runtime.py` caches the exact, independently hash-verified existing teacher and baseline arrays. The original scientific code remains unchanged; the computational amendment is recorded. A subsequent audit found repeated expression gating of V24 fallback coordinates. `v32_fallback_correction.py` replays the saved fits, preserving exact V24 outside supported outputs. Both original and corrected results remain available; only corrected results are interpreted.

## 2. V33: tabular foundation model

TabICLv2 is compared with ridge and CatBoost on IDENTICAL paired examples and features. This bounded pilot uses eight response programs, up to 512 route-balanced in-context examples, and two TabICL ensemble members. Own-output-excluded least-squares program coefficients avoid turning a missing own response into a supervised zero. The primary V24 correction strength is fixed at 0.5; strength 1 is descriptive only. An identity-shuffled TabICL is a negative control.

The public regression checkpoint is pinned to author repository revision `4dcd344ece2c00be9e831fdd35bed57b5ad83e19`. Its 114,324,594 bytes match LFS SHA256 `0db9cb538f114e79026bf08f45f41ad8dd7ad2de2aaca9a5ca8cd3bd9748ae7a`. BSD-3-Clause; ungated. The official loader uses `torch.load(weights_only=True)`. Dependencies are isolated inside `.v32_dependencies`, preserving the prior scientific runtime.

This is a small-budget local adaptation, not a reproduction of the authors' entire benchmark or a proof about all tabular model configurations. The V32 fallback correction is adopted before any V33 outcome is opened. The first adapter invocation stopped before scoring at a missing feature-scaling key; `run_v33_fixed.py` attaches the already fitted scaling values to the query adapter. Frozen initial code and both amendment records are retained.

## 3. V34: control-anchored conditional flow

Fit a residual flow directly on 1,024 RNA coordinates selected using allowed training sources. Condition on 32 STRING target coordinates and control-rank context features. Training matches unpaired, lane-centered control and treated populations through an OT assignment; it does not infer observed individual counterfactual pairs. Own-target expression is removed from coupling costs and endpoints/loss, preventing observed knockdown from serving as a shortcut.

The network has two 128-unit hidden layers. It trains for eight fixed epochs of 400 condition-balanced steps, using velocity, population-mean and direction losses. Select the epoch by eight-step GENERATIVE inference on source-only held-out target identities, rather than teacher-forced reconstruction. Prediction supplies no observed treated mean or source-treated posterior. Outside the modeled axis, retain V24. Pre-outcome code review also locked the untrained own-gene coordinate to its control value during every integration step; `v34_corrected_methods.py` preserves the initial implementation separately and prevents an unconstrained own-output head from changing later autoregressive inputs.

Training uses available raw H1/Orion/CD4 source cells, with the complete outer family excluded. K562/RPE1 pseudobulk means are not misrepresented as single-cell flow examples. V34 does not receive measured source-response tokens at inference, whereas V24/V32/V33 can use them. Its comparison with V24 changes both architecture and usable intervention information; failure does not isolate the flow model class. Missing STRING embeddings use a zero descriptor, a further limit of this pilot. Public mean comparisons use all previously frozen V25 cases. Fixed count diagnostics use the first four cases per context and one disjoint 400-cell control-null case per context; D1's 200-cell treatments are explicitly lower-power diagnostics. Compare unchanged V24, the learned mean through the same emitter, and a 25% flow-profile shape blend. Before any count outcomes, the shape comparison was amended to retain the learned-mean reference gene totals exactly, retain each cell library, and match its mean CPM probability within L1 0.0005. This separates variation from a simultaneous mean shift; the original frozen script remains available. The first integer projection failed the unchanged CPM tolerance in D1 (L1 0.0007056). `run_v34_projected_counts.py` applies the existing 16-depth-bin margin allocation uniformly to all cases; initial null outputs and the mechanical failure are retained. This also constrains gene-depth associations to the reference within those bins, a material limit of the shape test. Count metrics use the existing raw Wilcoxon/BH and panel-score adapter, without private host anchors.

## Reproduction

From Windows Command Prompt, in the project root. All Python commands below use:

```bat
set "SCI=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
```

The plans already exist and must not be overwritten. For a new clean experiment checkout, create plans with each script's `freeze` stage before running outcomes. Large arrays/weights remain untracked; their metadata/hashes are recorded.

```bat
"%SCI%" test_v32_transfer.py
"%SCI%" run_v32_cached.py
"%SCI%" v32_fallback_correction.py
"%SCI%" audit_v32_v34.py 32

"%SCI%" -m pip install --target .v32_dependencies --no-deps tabicl==2.2.0 einops==0.8.2 catboost==1.2.10
"%SCI%" collect_v33_checkpoint.py
"%SCI%" test_v33_tabular.py
"%SCI%" run_v33_corrected.py
"%SCI%" audit_v32_v34.py 33
"%SCI%" audit_v33_predictions.py

"%SCI%" test_v34_flow.py
"%SCI%" test_v34_integration.py
"%SCI%" run_v34_corrected.py train
"%SCI%" run_v34_corrected.py mean
"%SCI%" run_v34_projected_counts.py null
"%SCI%" run_v34_projected_counts.py treatment
"%SCI%" audit_v32_v34.py 34
```

To repeat the completed experiments without overwriting their frozen results, use the fresh-namespace runner (the original large input files are required):

```bat
"%SCI%" reproduce_v32_v34.py 32
"%SCI%" reproduce_v32_v34.py 33
"%SCI%" reproduce_v32_v34.py 34
```

It copies the identical plan into a timestamped `experiments/replays/` directory and runs the audited corrected implementations. This is an outcome-known numerical reproduction, not another independent validation. Completed results cannot be silently rerun into the same result namespace. The retained source/model hashes and row manifests support independent replay. Seed repeats are not independent biological studies.

## Sources

- [PerturbMap](https://arxiv.org/html/2607.28090v1): motivates paired response reuse, but its recipient-anchored fitting cannot directly serve a control-only unseen recipient. V32 learns a transferable context-conditioned map instead.
- [Tabular perturbation author code](https://github.com/royerlab/tfm-perturbation) and [TabICL implementation](https://github.com/soda-inria/tabicl): support testing this model class, not a private competition score claim.
- [GeneGeoFlow](https://arxiv.org/html/2608.06824v1): motivates control-anchored flow, OT coupling and direction learning. V34 does not reproduce its gene-geometry modules, CRISPRa benchmark, or claim identical generalization.