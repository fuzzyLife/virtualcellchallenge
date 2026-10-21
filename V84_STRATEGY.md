# V84: transfer response changes learned from paired experiments

V84 is a new, independently implemented response-transfer experiment. It uses measured responses to other interventions to learn how a response changes between contexts. It is not a reproduction of scDEBART, C3TL, or a claim about private leaderboard methods. See `V84_RESULTS.md` for the completed decision; this protocol preserves the hypothesis and fixed tests.

## Why this differs from recent submissions

V66 estimates a response from the same legacy K562, CD4, H1, HCT116 and HEK293T families used through much of this project. Many later versions change pooling, calibration, identity preservation or count realization rather than add independently measured response information. The latest published V82 package exactly reproduces V66's overall 0.1508821875498449.

The new Arc essential-screen tables provide no direct observations of our 300 targets. V83 rejected direct pooling on that coverage check. V84 instead uses 638 other interventions measured across K562, RPE1, HepG2 and Jurkat as paired examples for learning context changes. This tests whether useful cross-context response structure transfers to our different target panel.

## Data and provenance

- Arc dataset `arcinstitute/PIE_replogle_nadig_essential`, revision `91f86154aa1761a2dd0d885b979b9169f3b9983f`. Four full DE Parquet downloads were checked against publisher LFS SHA-256 values. Source README and inventories are retained in V83 analysis.
- Shared intervention eligibility: at least 30 perturbed cells and 400 reference cells in each of the four sources, determined from annotations before effect evaluation.
- Shared output axis: 6,087 mapped genes expressed at least 5 CPM in all four source reference profiles. All 300 competition genes are excluded as outputs of the new operator.
- Released reference cohorts differ by intervention. Their expression means are averaged for context routing; they are not counted as independent biological replicates.
- Released means are CP10k; multiplying by 100 gives CPM. Released log2 fold changes are used directly, then processed as specified below. No pseudo raw counts are inferred from these means.

## Fixed model

For each intervention, remove its own readout where present, apply the existing bounded knockdown-dose convention, clip log2 fold changes to +/-3, and subtract the row median on the shared output axis. The dose gain is `clip(0.9/max(knockdown,0.35),0.8,1.8)` when own-gene suppression is observed. This is a modeling assumption, not a measured efficiency for the competition contexts.

Split interventions deterministically using `sha256('v84-essential:'+target)%5`. The 513 fitting interventions determine all scales and parameters; 125 interventions remain held out. Scale each K562 output by the greater of its fitting-set standard deviation and 0.1.

The primary model uses a fixed rank-32 randomized SVD of standardized K562 responses. A zero-intercept ridge regression maps those scores to the standardized destination-minus-K562 response. The penalty is 0.1 times the mean diagonal of the score Gram matrix. The secondary model uses a separate zero-intercept slope per output gene, with a corresponding fixed relative penalty. Zero input gives zero correction.

For the negative control, refit destination responses after a cyclic derangement of the fitting target rows. This preserves source distributions and fitting machinery while destroying target pairing. Held targets never enter fitting. Also compare against copying the K562 response and against the destination fitting-set average response. No hyperparameters change after these comparisons.

To predict competition targets, use their existing measured K562 response on the shared axis. Transfer it through the fitted operators, then mix corrections using pairwise-valid Spearman correlations between query NTC expression and source reference expression. Prior weights are K562=1, RPE1=1, HepG2=0.5, Jurkat=0.5, multiplied by `exp(4*(rho-max(rho)))` and normalized; K562's correction is zero. Routing similarity is a heuristic, not evidence that response mechanisms match.

Add one quarter of the bounded +/-1 correction to the fold-specific V66 response. Own-gene responses stay zero because the unchanged count emitter handles own-gene suppression once. The existing 80% identity-margin guard is anchored at V66. Missing source coordinates abstain; absent query context correlations return the exact baseline.

## Development screen and separate donor test

Development uses all 655 previously reused public target/context pairs, excluding the recipient study from legacy response estimation. This is development, not blind validation. Require equal-study profile error / V66 <=0.995, worst-study ratio <=1.03, PDS loss <=0.005, and better profile error than the separately fitted shuffled model. No eligible arm means no selection. The prespecified primary `program25` is still donor-tested if development fails, with that failure retained.

The separate donor is Marson2025 D3 Rest, downloaded from CIPHER Zenodo record 21729034 and fully verified against MD5 `e8c066b72f12caea34d26985f49474ea`. It belongs to a previously used study, so this is a project-unseen donor, not an independent study or a final-test estimate.

Metadata selected all 14 panel targets with at least 400 supported cells: ADNP, AKAP1, MAP2K2, MARCKSL1, MYRF, PCGF1, PRDX3, SLC39A1, SUPT7L, USP28, ZC3HAV1, ZNF25, ZNF280C, ZNF514. The extraction has 18,463 cells and 18,533 aligned columns; 17,772 columns are measured in this source. Missing genes retain explicit measurement masks.

The 10,000 NTCs split within lane into 5,006 inference/template controls and 4,994 reference controls. Each target has 400 distinct actual cells and 400 distinct lane-matched NTC templates per seed. Seed offsets 0 and 620000 are technical repeats that may overlap cells. They are not separate biological experiments.

**Exclude every CD4 response and pooled donor summary.** Infer query expression using only the 5,006 NTCs. The CD4-excluded V24 baseline branch must not access the CD4 teacher. Independent poisoning checks alter excluded response values, uncertainties, doses, masks and context profiles and require bitwise-identical V24, V66 and V84 predictions. This donor comparator is the V66 algorithm with the whole CD4 family excluded, not the exact submitted ABC array.

The unchanged `residual10` count emitter uses matched templates, seeds and original NTC library sizes for baseline, primary and shuffled primary. Verify nonnegative integer counts, exact original row totals, eight independent tied-rank coordinates per population, BH correction, DE direction metrics, profiles and panel discrimination. Assess both seeds separately.

Count requirements retain the original thresholds: profile error must not increase; PDS may fall by at most 0.005, fidelity by 0.002, reach by 0.005, and Jaccard by 0.0005; available NMAE may increase by at most 1%. At least three endpoints must improve, and the genuine model must beat the shuffled model in profile error and either fidelity or Jaccard. Both development and donor requirements are needed for a release recommendation. Raw public scores are not converted into a fictitious hosted overall score.

## Conditional diagnosis if the primary fails

The oracle protocol was saved before reading donor scores. On six metadata/hash-selected targets, split the stored 400 actual cells into disjoint 200-cell fitting and scoring halves. Cross baseline versus measured-fitting-half means with old versus fitted population shapes. Within each mean pair, enforce exact gene totals, row totals and raw profile identity using the existing V29 projection. Record its additional CPM error contract instead of assuming it holds.

This asks whether supplying measured population shape is enough when the baseline mean stays fixed, and whether supplying the measured fitting-half mean makes a larger difference. It is not deployable: target treatment labels are available only to this diagnostic. It has noisy fitted means, 200-cell scoring power, only six targets, and depth/projection confounding when comparing across mean pairs. It must not train or select a new submission.

## Reproduce in a fresh analysis directory

Use the Python 3.12 scientific environment used by existing project analyses; `.competition_env` is the official CLI environment and is not interchangeable. Preserve existing frozen output directories rather than deleting them to rerun. Collector retries verify already downloaded files. Model fitting, extraction and inference deliberately refuse to overwrite some outputs.

```powershell
$sciencePython = 'C:/Users/animeshs/AppData/Local/Temp/vcc_v20_stack/env/Scripts/python.exe'
& $sciencePython -X utf8 collect_v84_sources.py
& $sciencePython -X utf8 v83_context_expansion.py sources
& $sciencePython -X utf8 collect_v25_holdout.py
& $sciencePython -X utf8 extract_v83_donor.py
& $sciencePython -X utf8 test_v84_operator.py
& $sciencePython -X utf8 v84_intervention_operator.py pairs
& $sciencePython -X utf8 v84_intervention_operator.py train
& $sciencePython -X utf8 audit_v84_operator.py
& $sciencePython -X utf8 run_v84_experiment.py prepare
& $sciencePython -X utf8 run_v84_experiment.py means
& $sciencePython -X utf8 run_v84_experiment.py donor
& $sciencePython -X utf8 audit_v84_exclusion.py
& $sciencePython -X utf8 finish_v84_donor.py
```

The committed inventories, plans and execution protocol are required inputs; download collectors do not invent or update their pinned revisions. `audit_v84_paper_code.py` separately checks reviewed author attention class definitions and architecture configuration. It does not load a checkpoint or implement the predictor above. Oracle execution, if its condition is met, uses `run_v84_oracle.py run` and then `assess`.

## Literature motivating the test

[scDEBART](https://proceedings.mlr.press/v306/sung26a.html) motivates basal-conditioned expression-change representations. Its reported fine-tuned benchmarks differ from the competition's controls-only unseen-context task. [Systema](https://www.nature.com/articles/s41587-025-02777-8) motivates the constant-response and shuffled-pairing checks. [C3TL](https://arxiv.org/html/2603.13051v1) motivates sharing intervention evidence across contexts; its context representation uses measured effects, so it cannot be assumed to work unchanged when only target-context NTCs are available. [Virtual Cells Need Context, Not Just Scale](https://pmc.ncbi.nlm.nih.gov/articles/PMC12919078/) motivates testing context diversity; it is a preprint, not a verified recipe for the competition.
