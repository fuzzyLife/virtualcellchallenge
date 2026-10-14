# V25: verifying cross-domain ideas before a new submission

The supplied bioCon1/bioCon2 posters motivated a new experiment. They are conceptual inputs, not instructions or evidence that their claims have been validated. Their file hashes are in analysis/poster_inputs.json. V24 remains the published reference:0.14299785521745675.

## What can and cannot be verified

| Poster idea | Evidence status | Consequence for this competition |
|---|---|---|
| Fisher information and natural gradients | Established mathematical tools; require a specified statistical model | Use regularization and a trust region to limit uncertain response corrections |
| Covariance and linear response | A useful predictive hypothesis under assumptions; not general causal identification | Anchor transfer to measured CRISPR effects and test it across held-out source families |
| Gaussian distribution transport | Well-defined mathematical maps under Gaussian assumptions | Test a bounded Bures covariance map against simpler controls |
| Exact universal neural/immune/epigenetic isomorphism | No verifiable primary paper found for the named MVGD framework or2040 experiment | Do not assume an identified cross-domain map or lossless sequence prediction |
| Equal dimension or Betti numbers imply a unique response map | Insufficient: topology does not identify gene correspondence, metric or intervention dynamics | Keep literal gene identities and verify alignment; do not infer response from visual similarity |
| Landauer/Carnot heat or entropy predictions from RNA snapshots | Required thermodynamic observables are absent | No calorimetry, heat, energetic-efficiency or physical-law claim |

[Martens' natural-gradient analysis](https://arxiv.org/abs/1412.1193) motivates curvature-aware optimization, damping and trust regions. It does not establish biological equivalence across substrates. [Transport information geometry](https://arxiv.org/abs/1803.06360) develops geometry of probability distributions; such a construction is not an empirical mapping from synaptic weights to immune receptor sequences.

[CIPHER](https://pmc.ncbi.nlm.nih.gov/articles/PMC12363937/) motivates examining co-fluctuation structure. Its [author forward implementation](https://github.com/GoyalLab/CIPHER/blob/main/cipher/forward.py) estimates a response amplitude using observed treatment shifts on selected training genes. That setting cannot be described as predicting an unknown query response amplitude from controls alone. V25 supplies source experimental responses and never reads query-treated cells during inference.

[TxPert](https://doi.org/10.1038/s41587-026-03113-4) combines basal-cell information with perturbation representations informed by knowledge graphs. It reinforces the separation of baseline state and intervention evidence. Its benchmarks do not establish performance on our leaderboard or a universal cross-domain law. Some publisher/PMC full-text requests were blocked; the primary abstracts, mathematical papers and author code were used, with this access limit retained.

## A direct falsification check

Consider a two-dimensional stochastic system:

```
dX = J X dt + sqrt(2) dW
J = [[-1, -omega], [omega, -1]]
```

For every omega, stationary covariance is I: J I + I J^T +2I=0. Both systems have the same stationary Gaussian distribution and the same Gaussian-mean Fisher metric. For a constant intervention force f=(-1,0), the mean shift is -J^-1 f:

| omega | Stationary covariance | Intervention mean shift |
|---|---|---|
|0|I|(-1,0)|
|2|I|(-0.2,-0.4)|

The identical snapshot geometry does not identify the intervention response. This is an explicit calculation, not a literature analogy. The rotation/nonequilibrium dynamics are unobserved. Three unit tests verify this counterexample, identity/Bures mapping, bounded singular values and disagreement abstention. The result limits our claim: V25 is a testable predictive approximation grounded in actual interventions, not a recovered physical law.

## Data audit uncovered a modality mismatch

The author's [Zenodo21729034 release](https://doi.org/10.5281/zenodo.21729034) includes curated fibroblast and RPE1 data. Metadata-only range inspection finds79 and34 overlapping competition targets respectively, but neither file has any overlapping target with200 or400 retained cells. Label namespaces are literal gene symbols, so this is not an Ensembl/symbol-mapping error. These are [Southard et al. CRISPR activation screens](https://www.nature.com/articles/s41588-025-02284-1), not knockdown labels. Activating a gene cannot generally be converted to knockdown by blindly negating its response. We therefore do not download8GB and silently add those responses to a CRISPRi teacher. Their controls and programs remain possible future representation sources, requiring separate cross-modality validation. RPE1/fibroblast compatibility would not prove that anonymous B/C have those literal identities.

A new Marson/CD4 D3 resting file has37 overlapping targets and14 groups with>=400 cells. It is selected by metadata before treatment scoring. The planned collector fetches bounded resumable ranges and checks the complete file against the author's MD5. The prepared raw auditor checks original cell IDs, sample lanes, labels, libraries and128 complete aligned source rows. In this run, all finalists failed development, so the incomplete download was stopped and these full-file/extraction/raw-row checks were not executed. Reserve targets are the first eight alphabetically meeting400-cell and lane-control support. No duplication manufactures power. The donor aliases cannot be matched confidently to the legacy CE donor IDs, so **the entire CD4 family is excluded from fitting on the reserve**, including all public CD4 summary effects and the V24 CD4 specialist. This is new raw donor data from an existing study, not a new independent biological study.

## New mechanism: intervention-aligned response geometry

V19 already tested covariance transport in16 control-variance PCA directions on2048 genes and failed release. Repeating that would not be a new idea. V25 instead fits64 directions to experimentally measured, training-only intervention mean shifts on up to4096 shared nonpanel output genes. It includes verified CD4 control geometries in source folds where CD4 is allowed. CD4 summaries and controls are condition/study-family matched; exact donor aliases are not established. This targets directions that actually respond to interventions, rather than directions with the largest spontaneous variation.

Source LFCs are corrected for measured source knockdown using the existing0.9 reference, with gain0.8..1.8, and clipped+/-3. Source means are mapped approximately into logCP10k response units with the derivative ln(2)*CP10k/(1+CP10k). This is a linearization, not exact reconstruction of a treated log-expression distribution. Missing source response coordinates are zero only inside the projection approximation if>=90% of the axis is available; poorer-coverage source targets abstain. They are not interpreted as biological zero responses. The300 target-panel genes are excluded from fitting the geometry.

Control expression is residualized for sequencing batch and log-depth using the existing verified V19 nuisance routine. Projected covariance is regularized:

```
C =0.75*empirical_covariance +0.25*diagonal +0.10*mean_diagonal*I
```

Two hypotheses transfer an observed source response z_s to query-control geometry:

```
force-response map:  C_query * inverse(C_source)
Gaussian Bures map:  C_source^-1/2 * (C_source^1/2*C_query*C_source^1/2)^1/2 * C_source^-1/2
```

The latter has a Gaussian optimal-transport interpretation before clipping. Operator singular values are clipped to2/3..1.5; clipping means we do not claim an exact covariance pushforward afterward. The force map assumes an approximately transferable effective force, a hypothesis the OU counterexample shows cannot be guaranteed by covariance.

Compute maps independently using two within-batch halves of the query's training controls. Opposing components abstain; magnitude disagreement shrinks corrections. Inverse regularized query covariance defines a Gaussian-mean Fisher proxy. A correction is limited to half the baseline response's norm in this proxy before its alpha and gene-coordinate bound. This is a mathematical stability limit, not a calibrated confidence interval. No query treatment information sets the amplitude or selects modes.

Add the bounded correction to the V24 response estimate and retain the same V16/residual10 raw-count emitter, target-own knockdown and exact template libraries. This phase tests response means, not an invented population-distribution improvement.

## Frozen tests and release policy

Finalists are force25, force50 and bures25. Controls: exact CV V24, diagonal25, alpha-matched shuffled25/shuffled50, amplitude25 and source_program25. The source-program arm tests direct smoothing/reconstruction without target covariance. Alpha is0.25 or0.5; every pre-gated gene correction is clipped to1LFC before alpha. The full plan and all numerical dependencies were hashed before outcome scoring.

Development reuses44 public cases:28 genuine400-cell cases in H1/HCT116/HEK293T, eight200-cell resting CD4 D1 cases and eight400-cell CD4 D2 stimulated48h cases. Those biology panels have been used before and are development only. Every evaluated family is excluded from fitted responses and source geometries; CD4 tests disable all CD4 mean teachers as well. The comparator is the V24 algorithm under the same source exclusions, not the exact published package. Complete HCT116 NMAE remains unavailable.

Apply the predeclared available-metric guards, require>=3 improvements including DE, and protect worst-family fidelity/reach. A finalist must beat its alpha-matched shuffle, the amplitude diagnostic and source-program smoothing on at least two fidelity/reach/Jaccard metrics. Source_program25 and amplitude25 are explicitly alpha0.25 diagnostics, so the force50 comparison is not perfectly amplitude-matched to these two controls. Only one winner selected on development can open the fresh donor reserve; no re-selection after reserve. Failure retains V24. A passing model still needs full production count/package/source-identity audits, and no hidden-score improvement is claimed before a new host result.

## Reproduce

Root scripts and small source/model/case/hash evidence are versioned; raw files, fitted arrays and count caches are ignored. Scientific runtime is the unchanged Python3.12.10 environment pinned by requirements-v23-science.txt, with two BLAS threads. Official packaging uses the competition environment separately. Existing plans/receipts refuse overwrite.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
set "VCC_SCI_PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
"%VCC_SCI_PY%" probe_v25_sources.py
"%VCC_SCI_PY%" audit_v25_source_labels.py
"%VCC_SCI_PY%" test_v25_geometry.py
"%VCC_SCI_PY%" prepare_v25_experiment.py
"%VCC_SCI_PY%" run_v25_experiment.py models
"%VCC_SCI_PY%" run_v25_experiment.py development
"%VCC_SCI_PY%" assess_v25_experiment.py development
"%VCC_SCI_PY%" collect_v25_holdout.py
"%VCC_SCI_PY%" extract_v25_holdout.py
"%VCC_SCI_PY%" audit_v25_cd4_inputs.py
```

Only a passing locked finalist opens confirmation. Actual results and the final release decision are recorded separately in V25_RESULTS.md and V25_SUBMISSION.md. No failed historical threshold is relaxed.

Before opening development outcome tables, release_guard_amendment.json added conservative diagnostics: force50 must also beat separately emitted alpha0.5 amplitude/program controls, and bures25 must beat its own Bures-coordinate shuffle on two DE metrics. These can only reject the selected candidate; no fallback or re-selection is permitted. No model, cell, seed or original threshold changed.
