# V59: state composition and conditional response

This experiment tests the user's requested new release. V24 remains the measured comparator until a new candidate passes. No previous code, prediction package or frozen model is modified. Public development biology has already been examined; generation seeds are not biological replicates.

We distinguish changes in the relative frequency of RNA neighborhoods from expression changes within those neighborhoods. V20/V21 used a response-driven occupancy tilt; V59 instead estimates abundance and conditional responses directly from public experimental populations. Four states are shared across source studies through a training-only eight-coordinate PCA. The 512 input genes exclude all 300 panel targets and are selected for low perturbation response in allowed training sources. NTC cells alone determine each lane's normalization. Treated RNA can still change a cell's assigned state; this is a predictive decomposition, not an identified survival or lineage mechanism.

The complete recipient study is excluded from stable-gene selection, PCA, clustering and donor responses. HCT116/HEK293T are one Orion study; resting/stimulated CD4 share one total source vote. Query controls calibrate their own lane transforms. Missing genes/states abstain. Matching source interventions are valid cross-context evidence; this is not strict unseen-target evaluation. D3 remains unopened.

First reconstruct an independent treatment split from the other split's estimated density ratio and conditional response. Compare composition-only, within-only, joint, direct mean and state-permuted controls. These oracle diagnostics cannot enter submission predictions. All raw row IDs come from unchanged V52 treatment/control views. Oracle results average lanes and views within target before averaging targets.

The deployable primary is fixed joint50: a 50% bounded correction toward the new composition/within-state response, retaining V24 outside supported coordinates. Diagnostic arms are composition50, pooled50 without query abundance changes, and a within-source target-deranged joint50. All use saved V58 comparator arrays (H1/Orion comparators additionally exclude CD4). Thus the clean held-CD4 comparator is unchanged, while H1/Orion results concern this stricter baseline.

Count tests use 400 generated cells and up to 400 actual distinct treatments with independent real NTC references. Control pools/references are common across targets per source/seed for interpretable panel discrimination. All five mean arms use identical templates and the unchanged residual10 emitter. A sixth diagnostic uses a half-and-half mixture of V24-edited and learned-state-weighted real NTC cells, with stochastic balanced edits and exact libraries; it does not enforce V24 aggregate means. Sampling is without replacement. Shape changes are evaluated after integer generation.

The release gates are frozen in analysis/plan.json before evaluation. Response prediction must beat V24 by at least 1% in equal-study MSE, avoid a >2% worse study, and beat shuffled and pooled controls in each study. Counts must improve expression and at least two DE metrics, satisfy fidelity/Jaccard guards and beat shuffle. A diagnostic cannot be promoted after seeing outcomes. A failed candidate is documented and cannot be called an improved submission. No automatic leaderboard upload or Git push occurs.

Primary methodological references: [Milo](https://www.nature.com/articles/s41587-021-01033-z) distinguishes differential abundance on cell neighborhoods; [muscat](https://pubmed.ncbi.nlm.nih.gov/33257685/) studies subpopulation-specific differential state. Neither establishes VCC performance for this proposal.

## Reproduce (Windows CMD)

Use the existing pinned scientific environment; the competition environment is reserved for official packaging.

```bat
cd /d L:\promec\Animesh\Download\Virtual
set "SCIENCE_PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
"%SCIENCE_PY%" -X utf8 test_v59_state_composition.py
"%SCIENCE_PY%" -X utf8 run_v59_experiment.py response
"%SCIENCE_PY%" -X utf8 run_v59_experiment.py counts
```

Plans refuse overwrite; completed immutable stage artifacts are verified and resumed. Root scripts and compact receipts are versioned. Large matrices are regenerable local caches under the system temporary directory. The final result report records what actually completed and whether a submission was justified.
