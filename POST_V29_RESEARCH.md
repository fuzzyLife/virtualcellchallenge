# What remains untested after V29: leaderboard and public-method review

Review date: 2026-10-15. This is research and code inspection, not a new performance result. **V24 remains the best published project entry, overall 0.142997855.** No model was trained, reserve opened, package built or submission uploaded in this review.

The strongest next experiment is a broader, biologically informed **mean-response predictor** with consistent fitting/query inputs. Run a compact STRING-bilinear baseline before a large new generative network. Later couple an improved mean model to population transport. We have already tried neural mean heads; the proposal is a specific new representation, training-data scope and supervision, not our first intervention model.

## 1. What the leaderboard actually tells us

The supplied screenshot shows rank524 among1537 teams, overall0.143. Its rows are the teams' latest entries; historical ranks are not a model comparison. Our exact V24 receipt supplies the six scaled and raw scores. Competitor values are manually transcribed from the screenshot, rounded, and not independently retrieved from the API. Their anchor IDs are not visible. The public leaderboard HTML did not expose its table to the web reader.

The deterministic comparison is [metric_gaps.md](analysis/post_v29_research/metric_gaps.md); numerical provenance is [gap_summary.json](analysis/post_v29_research/gap_summary.json).

| Scaled component, higher is better | V24 | Illumina AI / iCell | X-team / X-89 |
|---|---:|---:|---:|
| Discrimination |0.6442|0.874|0.776|
| Expression accuracy |0.0476|0.787|0.458|
| DE fold-change accuracy |0.0823|0.419|0.283|
| Direction fidelity |-0.0158|0.153|0.045|
| Direction reach |0.0848|0.444|0.362|
| Significance overlap |0.0149|0.088|0.088|

Expression and fold-change accuracy account arithmetically for about56% of the rounded component gap to iCell. The third-place entry has overlap0.026, fairly close to our0.0149, while expression accuracy is0.455. Recovering more significant genes alone is therefore an incomplete target.

Raw errors tell the same story: V24 expression error is0.960415 and fold-change error0.950686, versus iCell's screenshot0.240 and0.748. Scaled NMAE0.0823 is a higher-is-better improvement score, **not an8.2% raw error**. This comparison is descriptive: it does not establish which architecture, data or causal mechanism produced another team's score.

[Official scoring specification](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics.md) documents the six equally weighted, reference-anchor-scaled components. The deployed submission contract and current FAQ take priority over older examples in repository documentation. Unanchored public-development metrics cannot be averaged into an estimated official score.

## 2. Public 2026 methods found, and the evidence limits

Searches covered the named leaders, their exact model names, public GitHub repositories and participant write-ups. **No verifiable architecture/training recipe tied to the screenshot's iCell, X-89, m0983n or GG-21 was found.** This is a search limitation, not proof that none exists. Model names are not reliable architecture descriptions; iCell is not established to be an unrelated similarly named published model. The [official challenge](https://virtualcellchallenge.org/) explicitly allows public or private data and does not require open code or weights. We cannot assume all teams have identical data, or assert that any particular private dataset was used.

- [Sidechain's participant log](https://saberhq.com/sidechain/models/) reports SER-17 at0.1657, submitted October8. Its disclosed design pools measured intervention effects, uses STRING neighbors and post-pooling shrinkage, and controls count emission. A real-control emission variant reportedly fell to-0.0039; selectively narrowing predicted variation recovered much of the loss. These are author-reported results, not independent reproduction or a top-team recipe. Many underlying ideas overlap V10-V16; this is not evidence that simply copying its emitter reaches0.46.
- [mmyen's 2026 repository](https://github.com/mmyen/virtual-cell-challenge-2026) implements a principal-component bilinear ridge predictor pooling Replogle/Nadig effects. Its README lists SLIM as **not implemented** and does not establish a competitive leaderboard score. The inspected implementation assumes particular preprocessed input matrices: its preprocessing should not be copied onto our raw-count files without auditing units and controls.
- [pphouse's 2026 repository](https://github.com/pphouse/virtual_cell) describes STRING-based response transfer and low-rank prediction across H1/K562/RPE1. Its documented completed early submissions have negative scores; a later entry is listed as scoring. It provides useful diagnostics, not verified evidence of a high-ranked solution. Its stated overlap pertains to its own selected library, not ours.
- [CaiZongyuan's 2026 repository](https://github.com/CaiZongyuan/virtual-cell-2026) reports that its first STATE fine-tune was worse than NTC and pauses submissions for local diagnosis. A foundation-model name alone does not imply good transfer.
- [Kaipeng Zheng's AtlasShift repository](https://github.com/kaipengm2/Virtual-Cell-Challenge-2026) reports0.1545618019 and rank82 **on September10**. Its historical rank is not a current top100 claim. We already reviewed and incorporated related atlas-transfer ideas before V5.

Primary code revisions and six file hashes are retained in [sources.json](analysis/post_v29_research/sources.json). External code was read only; no imports, installs or execution of that code occurred. Selected fitting methods and documentation were inspected, not every file in these repositories. Search-index summaries, self-reported results, and our reproduced experiments are separate evidence categories.

## 3. Where our experiments leave a real gap

The important local sources inspected are `train_v7_transition.py`, `train_v26_decoder.py`, `v20_methods.py`, V14/V18/V27/V29 results and the two retrospectives. Their hashes are archived in the gap summary.

- **We have trained intervention mean predictors.** V7 used target identity,16 Tahoe components and context ranks, with linear/64-unit/128-64-unit heads. Its inputs and source-response cache focus on the300 challenge targets. V14 added bounded context interactions; V18 evaluated pretrained PIE; V20 evaluated Stack-derived response experts. These failed or produced small gains. Recommending just another MLP or ensemble is repetition.
- **V26/V27 did not replace the mean predictor.** Training supplies observed source fitting-group means as decoder bases and descriptors derived from observed source effects. Query means remain V24 and query descriptors come from its predictions. Consequently, better reconstruction ELBO cannot demonstrate better query mean prediction. Source/query condition distributions also differ. This motivates consistent input conditioning, but does not establish a unique explanation for failure.
- **Distribution repair worked numerically and failed biologically.** V29 retained V24 means and achieved0-1 null DE calls across15 tests, yet H1 fidelity fell0.647 to0.102 and overlap0.057 to0.025. A better null is necessary evidence about calibration, not sufficient evidence about treatments.
- **Observed response estimates are weak.** On four already-spent H1 targets, predicted absolute effect sums on true DE coordinates are19-30% of observed sums, with sign agreement61-75%. This is selected-on-true-DE, development-only evidence, not a universal competition calibration. Oracle scaling helped some targets and worsened ACLY: multiplying every effect is unsupported.
- **Existing ensembles share errors.** V15 experts had residual correlations around0.99-1. Voting cannot introduce intervention information absent from every member.

## 4. Highest-priority untested experiments

### A. Fit a full-public-panel STRING-bilinear mean expert

[SLIM's August2026 preprint](https://www.biorxiv.org/content/10.64898/2026.08.07.743481v1) and [author code](https://github.com/RasmussenLab/SLIM) use STRING perturbation embeddings, an output response basis and a ridge-fitted bilinear map. The reported evaluations are predominantly within-dataset perturbation benchmarks, not a demonstrated solution to our unseen-context raw-count task. Selected author `BilinearModel.fit` code was inspected at revision5a7e9ade5d0a6b6331e6dbc81181450605047bcc.

**New here:** fit this explicit compact baseline on all suitable public perturbation targets, including auxiliary non-panel targets; do not restrict the fitting dictionary to300 target IDs. Preserve reliable direct-response evidence as an expert, and learn a bounded correction to it. Compare against a purely shared bilinear expert before adding low-rank interactions with NTC context features. Our previous low-rank programs, Tahoe heads and STRING-informed pretrained PIE inference are related but do not constitute this fit.

Why auxiliary data can matter: Replogle RPE1 was already audited as having zero direct target overlap. It cannot be a direct donor for these300 genes, but responses to other targets can train a structured target-to-response relationship. The value of this indirect supervision is **untested**, not guaranteed. Likewise HepG2/Jurkat evidence used by a pretrained checkpoint is not equivalent to fitting our own full-axis multi-context teacher on its audited data.

Audit the exact available targets, measured output genes, duplicate cells/studies, raw versus normalized matrices and guide metadata first. Build masked source-specific response summaries with training-only bases. Exclude each example's own target from response fitting; trans effects on other panel genes remain legitimate outputs. Keep missing measured coordinates missing. Fit pooled-profile effects and mean-cell-composition effects explicitly rather than treating them as interchangeable.

Start with a small regularization/basis grid. Required controls: V24, a shared mean response, an amplitude-matched baseline, random/shuffled target embeddings and the same model without context features. Novelty must survive those controls. First use the unchanged V24 emitter to isolate mean gains; then compare V29's calibrated emission in a separately declared mean-by-shape factorial.

### B. Re-estimate source effects with guide/assay-aware supervision

A source effect is currently often pooled before a model can distinguish guide efficiency, assay, batch and context response. V22 learned continuous fingerprint proxies but explicitly did not fit a separate guide-efficiency parameter. Its limited state interactions and V17 scalar mixture are not this experiment.

Fit a hierarchical source-response model using multiple guides/constructs where actually available, assay metadata and measured own-target suppression with uncertainty. Use non-target downstream genes to estimate response programs, and compare ordinary pseudobulk against an author RNA-fingerprinting teacher where feasible. The [RNA-fingerprinting author implementation](https://github.com/satijalab/rna-fingerprinting) had been collected earlier but not fitted as that teacher.

Missing guide metadata must trigger abstention or a declared simpler model, not fabricated guide labels. The query's hidden efficacy is unobserved; use a source-learned, assay-aware prior and sensitivity analysis rather than claiming it can be recovered from NTCs. Source target suppression may supervise training even though the own gene is excluded from leaderboard scores. Never use observed held-query treatments to fit efficacy or amplitudes.

### C. Train conditional population transport with a learnable mean

[Conditional Monge Gap](https://arxiv.org/abs/2504.08328) learns conditional transport from unpaired populations; [CellFlow's author tutorial](https://cellflow.readthedocs.io/en/latest/notebooks/100_pbmc.html) distinguishes deterministic flow matching and stochastic GENOT. These were researched but never trained here. V25's covariance/Bures rule, V29's integer margin repair and the V26 logistic-normal decoder are different methods.

Train from public control populations to public treated populations using target embeddings, control-derived context/state and measured source assay features. The model must predict the treated mean rather than receive it as an input. If a predicted-response descriptor is used, generate training descriptors out of family, so fitting and inference see similarly estimated inputs. Compare with an ordinary residual mean head and shuffled condition labels. OT pairs are modelling couplings, not observed same-cell counterfactual pairs. Published drug/cytokine generalization does not prove CRISPRi unseen-lineage transfer.

This is a larger, higher-risk experiment. Use a modest latent pilot and a tested count adapter; train it only after the mean-only baseline establishes that it contributes new response information.

### D. Full pretrained transition fine-tuning, rather than frozen-feature blending

A genetic-perturbation [STATE transition model](https://github.com/ArcInstitute/state), or an eligible public [TxPert graph model](https://github.com/valence-labs/TxPert), is a separate future experiment from our small V7 heads and frozen Stack/PIE blends. We have not run a complete target/context transition fine-tune with our raw-count release contract. This is partially explored, not an untouched foundation-model hypothesis.

TxPert's author README distinguishes public graph configurations from stronger proprietary-graph variants. Its processed5,000-gene/log-normalized output is not our18,533-gene raw-count submission. Audit vocabulary, license, preprocessing and checkpoint training families before a run; do not claim a public checkpoint reproduces private-model paper results. The public failed STATE example above is a reason to demand matched baselines, not abandon or trust the whole class.

## 5. Validation order and decision

1. Freeze the mean-model inputs, data units, training exclusions, candidate grid and negative controls. All H1 development cases discussed above are already spent. Repeated seeds improve sampling diagnostics but are not new biological contexts.
2. Use nested whole-study/family holdouts for model selection. Fit graph/output/context preprocessing on allowed training data only. Include both known-target/new-context transfer and a separate withheld-target stress test. Do not prohibit public evidence for a known target in the primary competition-like setting.
3. Measure pooled expression error and per-cell mean-composition fold-change error before generating360,000 cells. Require consistent benefit over V24 and matched shuffles, with family-level reporting. Report missing/low-power metrics rather than inserting zeros.
4. Test only mean finalists through400-cell null and treatment comparisons with the unchanged V24 emitter. A second, declared factorial can test mean/calibrated-shape interaction. Preserve discrimination/fidelity/overlap guards; no local raw-score average is presented as a leaderboard forecast.
5. Only after development criteria pass, complete and verify the untouched donor reserve. Exclude all related CD4 evidence from fitting during its evaluation. A new donor is additional confirmation, not an independent-study guarantee. No reserve is opened in this research turn.

**Recommended first action: A, then B.** They target the largest observed metric gap and add different supervision before adding generation complexity. Further majority voting, universal amplification, broad variance restoration and control-only geometry have already received substantial tests and are lower priority.

## Reproduction

From Windows Command Prompt in the project root:

```bat
.\.competition_env\Scripts\python.exe review_post_v29_leaderboard.py
.\.competition_env\Scripts\python.exe review_post_v29_leaderboard.py --check
```

The script uses only the Python standard library. It verifies V24's six-component arithmetic, checks screenshot rounding consistency, records local evidence hashes, and deterministically writes the metric table and JSON gap analysis. `sources.json` records pinned public code URLs and decoded UTF-8 hashes. The screenshot itself remains at its user-supplied location and is identified by hash; it is not copied into Git. This verifies traceability and arithmetic, not any proposed model's performance.
