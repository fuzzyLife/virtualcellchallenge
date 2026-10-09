"""Document post-V11 research, measured diagnostics and frozen future experiments."""
import json,sys,importlib.metadata
from pathlib import Path
import numpy as np,pandas as pd
import build_submission as build
from record_v10_result import table
ROOT=build.ROOT;OUT=ROOT/'analysis/post_v11_research'

def main():
 neutral=pd.read_csv(OUT/'neutral_summary.csv');n=neutral[neutral.reference_cells==2000][['method','mean_null_DE_calls','max_null_DE_calls','mean_variance_ratio','mean_detection_MAE']]
 transfer=pd.read_csv(OUT/'transfer_family_summary.csv');coverage=pd.read_csv(OUT/'source_target_coverage.csv')
 assessment=json.loads((OUT/'viperturb_assessment.json').read_text());plan=json.loads((OUT/'future_experiment_plan.json').read_text());literature=json.loads((OUT/'literature_review.json').read_text())
 commit=json.loads((OUT/'data/fingerprint_commit.json').read_text())['sha'];tm=json.loads((OUT/'transfer_audit_manifest.json').read_text())
 neutral_records=pd.read_csv(OUT/'neutral_scores.csv');assert len(neutral_records)==48 and neutral_records.exact_original_depths.all() and neutral_records.positive_integer_counts.all()
 assert set(neutral_records.context)==set('ABC') and len(neutral_records.seed.unique())==2
 baseline=json.loads((ROOT/'experiments/v11_consensus_coverage/analysis/reference_model.json').read_text())
 for name,h in baseline['dependencies'].items():assert build.sha256(ROOT/name)==h,name
 frozen=json.loads((ROOT/'experiments/v11_consensus_coverage/analysis/publication_manifest.json').read_text())
 for name,h in frozen['generation_snapshot_sha256'].items():assert build.sha256(ROOT/'experiments/v11_consensus_coverage/analysis'/name)==h,name
 sources='\n'.join(f"| [{s['title']}]({s['url']}) | {s['project_implication']} |" for s in literature['sources'])
 future='\n'.join(f"| {t['id']} | {t['hypothesis']} | {t['design']} | {t['negative_control']} |" for t in plan['tests'])
 confirmation=pd.DataFrame([{'context':c['context'],'targets':', '.join(c['targets']),'available_cells':', '.join(str(c['available_treated_cells'][t]) for t in c['targets'])} for c in plan['new_count_confirmation']])
 comp=pd.read_csv(OUT/'on_target_composition.csv');top=comp.loc[comp.other_gene_log2FC_at_kd90.idxmax()]
 text=f"""# Post-V11 research and tests

Reviewed and tested on 2026-10-09. V11 remains the published baseline: overall 0.137003974, entry `bmIAGveuEvFiYQTvz0Ma`. This work collects small public inputs, diagnoses remaining failures and freezes future tests. It creates no upload candidate and changes no V11 numerical code, model, package or frozen evidence.

## Main findings

The highest-priority missing opportunities are an independently profiled, probe-based genetic source; gene-specific count-distribution calibration; and calibrated context-transfer uncertainty. More elaborate networks are a later option, after those comparisons. The new neutral tests support the dispersion hypothesis but do not establish improved perturbation or leaderboard performance.

1. The original VIPerturb-seq target manifest contains **{assessment['manifest_exact_target_overlap']}/300 targets**. Only **{assessment['qc_filtered_panel_targets']}** are included in its filtered object, the median available count is **{assessment['median_available_target_cells']:.0f} cells per target**, and none has 400 cells. This may improve source/platform evidence, not solve validation power.
2. Across six A/B/C control splits, gene-specific gamma-Poisson emission reduces mean null DE calls against 2,000 independent references from **{float(n[n.method=='v11_noise65'].mean_null_DE_calls.iloc[0]):.1f} to {float(n[n.method=='gene_gamma_poisson'].mean_null_DE_calls.iloc[0]):.1f}**. Average total CPM variance moves from 0.233 times real-control variance to 1.051. It still makes up to 291 null calls in one split, so it is not fully calibrated.
3. A 25% empirical-control template blend worsens the null bias despite matching population moments. This tested blend should not replace V11. The participant Sidechain likewise reports that another real-control emission change failed on the board after local improvements; the implementations differ and this is not an independent reproduction of its experiment.
4. Whole-family summary diagnostics find weak target-specific transfer and poor reliability in the strongest heuristic-confidence bins. The current scalar pooling uncertainty is not a calibrated directional posterior. Source masks also leave a median of **{coverage.outputs_no_source.median():.0f}/18,532 non-target outputs** without a usable response per target. These masks include abundance/QC restrictions, not just unmeasured assay genes.

## Scientific papers and public evidence

The source record in `analysis/post_v11_research/literature_review.json` states exactly which sections or primary excerpts were reviewed. Full publisher pages for some Nature papers were inaccessible; the review uses primary abstracts/method excerpts and author code, and does not claim every supplement was read. VIPerturb-seq and the pinned RNA-fingerprinting source describe preprint work. Arc's experiment article and Sidechain's participant log are listed separately from scientific papers. We cannot infer undisclosed methods of higher-ranked teams.

| Primary source | What it changes in our test design |
| --- | --- |
{sources}

## Collected data and method code

The root `collect_post_v11_data.py` performs bounded unauthenticated downloads and verifies hashes when reusing them. The original Zenodo manifest additionally matches the publisher MD5. URLs, retrieval times, sizes, hashes and publisher checksums are in `analysis/post_v11_research/data_manifest.json`.

- **VIPerturb-seq:** original record `10.5281/zenodo.18460279`, manifest and target-level coverage/cell/QC table. The genome-wide experiment is K562; K562, HEK293 and HAP1 also appear in a separate small pilot. The genome-wide source must be grouped with other K562 data when holding out a lineage. The three full bins total about 10.2 GB compressed; their shared NTCs must be deduplicated. We have not downloaded the expression objects, verified their gene/count axes or extracted responses. The supplied Seurat RDS format needs an R/Seurat extraction route; no R runtime is available here. The QC-filtered object alone cannot supply all 300 targets. Do not substitute the VIM-enriched screen for the unbiased genome-wide screen.
- **Reactome:** exact downloaded GMT ZIP, derived pathway membership and target annotation coverage. The contemporaneous version endpoint reports **97**; the ZIP is pinned by bytes/hash rather than assuming its mutable URL guarantees a release. There are **2,868 human pathways**, and **184/300 targets** have annotations. Membership is a prior, not an intervention label. Data are CC0 under the [Reactome license](https://reactome.org/license); attribution is retained here.
- **RNA fingerprinting:** author README, MIT license and `R/fit.R`/`R/utilities.R` pinned at commit `{commit}`. The inspected implementation estimates 30 control PCs, removes baseline structure, then estimates a rank-one perturbation factor and filters it. It has not been installed or run. This is materially different from our existing control reconstruction networks and pooled pseudobulk effects. Mapping an observed treated query to a dictionary cannot directly predict an unobserved challenge response from controls alone.
- **Existing public sources:** original Replogle author inventory was refreshed; the V7 cached source response arrays are hash-verified and reused without changes. No newly verified independent lineage has been added to training.

See `viperturb_target_coverage.csv`, `viperturb_inventory.csv`, `reactome_pathways.csv`, `reactome_membership.csv` and `target_annotation_coverage.csv` in the research directory. VIPerturb data are CC BY 4.0 according to the original record; cite Bradu and colleagues and the dataset DOI when using them. Source code is retained with its MIT license. The data collection contains no challenge treated labels.

## Tests executed: distribution and null bias

`benchmark_post_v11_distributions.py` freezes `neutral_test_plan.json` before execution. Every guide is split into disjoint fit, template and reference populations (50/25/25%). All candidates use the same 400 template depths; references contain 400 or 2,000 distinct cells. No own-target knockdown or response is applied. There are four methods, three contexts, two seeds and two reference sizes: **48 comparisons**. All generated counts are positive stored integers and preserve original library depths exactly.

| Method | Description |
| --- | --- |
| real_resample | Independent real control cells; a benchmark, not moment-fitted |
| v11_noise65 | Exact V11 zero-response emitter and integer dual-moment fitting |
| gene_gamma_poisson | Per-gene latent overdispersion estimated from fit controls after counting-noise correction; independent gamma-Poisson template followed by the same moment solver |
| empirical_blend25 | 25% real-cell composition and 75% compact V11 template, followed by the same moment solver |

The two pilots keep requested mean and pooled compositions fixed. Their initialization changes; this is not a trained scVI or CellOT implementation. Independent gamma-Poisson draws do not retain cross-gene biological covariance. The empirical blend's weak variance and high null DE show that adding real-cell structure before moment fitting does not guarantee preservation afterward.

{table(n)}
False/null calls here are genes declared DE by the existing asymptotic Wilcoxon/BH implementation with a shared reference CPM filter greater than 5. Detection MAE is the average absolute gene-detection-fraction error on reference-eligible genes. Total variance matching alone is insufficient: B still has substantial null bias under gamma-Poisson. The per-gene diagnostic CSVs show where variance and detection remain wrong. There is no perturbation-level or official anchored score for these new pilots.

## Tests executed: transfer, uncertainty and systematic signal

`audit_post_v11_transfer.py` evaluates **{tm['rows']} public source contrasts** from five families/seven contexts, excluding the entire evaluated family from pooled responses and shrinkage. It reproduces the ungated V11 consensus-20 response equation; it does not apply the production context-abundance gate or generate counts. Unsupported contrasts abstain (**{tm['unsupported_contrasts']}** have no evaluable downstream outputs).

Matched-target, generic training-centroid and shuffled-target cosine comparisons use the same common predictor-supported output gene mask in each context and exclude all 300 panel outputs. Two target-label permutations test identity signal. Separate descriptive centering subtracts evaluation-panel centroids solely while calculating a residual metric; no held response enters fitting. This centered metric complements the official common-control metric and is not a replacement scorer.

{table(transfer[['family','lfc_mse_ratio','de_sign_accuracy','heuristic_95_interval_coverage','transfer_cosine','residual_cosine_diagnostic','shuffled_target_cosine','generic_train_centroid_cosine']])}
These are retrospective, unanchored source-summary diagnostics. Normal/SE source q-values define the reference-DE set; they are not generated-count Wilcoxon calls. The MSE ratios above 1 show that this ungated estimator can be worse than zero effects in fold-change MSE even where some directions transfer. This is not a contradiction of V11's published score: the representation, gating and metric are different. A generic shared effect can have a larger cosine while containing no target-discrimination information.

`uncertainty_bins.csv` separates support, agreement and effect/uncertainty bins. The heuristic uncertainty is pooled transfer variance plus held-source measurement error, not an empirical-Bayes posterior variance. Its normal sign-error proxy and observed sign errors among reference DE genes should therefore be compared descriptively, not treated as a validated false-sign posterior or unconditional calibration curve. The very strong bin has many more errors and much poorer interval coverage than its proxy suggests. Future calibration must be fitted on inner training folds and confirmed independently.

## Small hypothesis ruled down: target knockdown composition

`audit_post_v11_composition.py` tests the exact mean-composition algebra for changing only the target: other-gene log2 fold change is `-log2(1 - knockdown * target_fraction)`. The largest 90%-knockdown effect is **{top.other_gene_log2FC_at_kd90:.6f}**, for `{top.target_gene}` in `{top.context}`; changing 80% to 90% changes it by only **{top.other_gene_log2FC_kd90_minus_kd80:.6f}**. This is a real indirect effect even though the own target is excluded from scoring, but its small size makes it a lower priority than the measured distribution and source-transfer gaps. This algebra is not a test of absolute RNA abundance or statistical significance.

## Concrete next experiments, in order

The frozen `future_experiment_plan.json` contains the full designs, negative controls and screening limits. All comparisons retain exact V11 as a baseline. Source additions, response estimation and count emission must be isolated before combining them.

| Test | Hypothesis | Design | Negative control |
| --- | --- | --- | --- |
{future}

There are 12 reserved, globally disjoint, previously unused count-test targets:

{table(confirmation)}
Selection used target names and available cell counts in existing control-defined technical strata, not response values. All previous V7-V11 count-test targets are excluded. Public source summaries have already been inspected, so this is not pristine biological blind validation. H1 remains reused development data. Preserve genuine observed counts (often 150-200), and never duplicate cells to claim 400 independent observations. Use both sampled-400 and full independent references, paired seeds and equal family/technical-stratum weighting.

The initial count-stage is limited to four candidates including V11. Freeze the development winner before opening reserved confirmation. Expression error must improve on both stages, while fixed guards protect fidelity, reach, fold-change accuracy and discrimination; inspect the worst family as well as the mean. A failed confirmation must remain a failure. These guards are screening criteria, not a forecast of the official anchored score. A neutral-only success does not authorize claiming an improved submission.

## Reproduce from Windows Command Prompt

```bat
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe collect_post_v11_data.py
.\\.competition_env\\Scripts\\python.exe test_post_v11_research.py
.\\.competition_env\\Scripts\\python.exe audit_post_v11_transfer.py
.\\.competition_env\\Scripts\\python.exe benchmark_post_v11_distributions.py
.\\.competition_env\\Scripts\\python.exe audit_post_v11_composition.py
.\\.competition_env\\Scripts\\python.exe document_post_v11_research.py
```

`plan_post_v11_experiments.py` was run once before future modeling. Its frozen JSON is tracked; rerunning it deliberately refuses to overwrite the plan. Neutral-test checkpoints verify plan and per-gene-file hashes before reuse. Data downloads verify cached bytes; changed upstream inputs require an explicit new snapshot rather than silently replacing this one. New code uses the existing competition environment and leaves its dependencies unchanged.

`test_post_v11_research.py` passes held-family leakage, masked-value abstention, no-source abstention, cosine missingness, count/depth-axis and valid-template checks. Documentation requires all 48 neutral records, exact-depth/count invariants and unchanged hashes for V11's numerical dependencies and 17 frozen generation snapshots. `research_manifest.json` records the complete input/output/code provenance. No submission was generated or uploaded by this work.
"""
 (ROOT/'POST_V11_RESEARCH.md').write_text(text,encoding='utf-8')
 readme=ROOT/'README.md';s=readme.read_text(encoding='utf-8');marker='<!-- post-v11-research:begin -->';end='<!-- post-v11-research:end -->'
 if marker in s:s=s[:s.index(marker)]+s[s.index(end)+len(end):].lstrip('\n')
 block=marker+'\n\n## Post-V11 research and tests\n\nReviewed ten relevant scientific papers/methods and separately checked official data and participant reports. Collected the VIPerturb manifest (300/300 targets; 124 QC-filtered; median 47 cells), pinned Reactome annotations and RNA-fingerprinting code, and ran 48 control-distribution comparisons plus five-family transfer diagnostics. A gene-specific gamma-Poisson pilot reduces average null DE calls from 1,711 to 84 against 2,000 references, but retains context-specific bias and has no perturbation/leaderboard validation. Twelve unused count-confirmation targets are frozen for the next experiments. V11 remains the best submitted model. See [POST_V11_RESEARCH.md](POST_V11_RESEARCH.md) for papers, findings, limitations and reproduction.\n\n'+end+'\n\n'
 s=s.replace('## Completed work',block+'## Completed work',1);readme.write_text(s,encoding='utf-8')
 files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='research_manifest.json')
 scripts=['collect_post_v11_data.py','audit_post_v11_transfer.py','benchmark_post_v11_distributions.py','audit_post_v11_composition.py','plan_post_v11_experiments.py','test_post_v11_research.py','document_post_v11_research.py']
 build.write_json(OUT/'research_manifest.json',{'report_sha256':build.sha256(ROOT/'POST_V11_RESEARCH.md'),
  'script_sha256':{n:build.sha256(ROOT/n) for n in scripts},'artifact_sha256':{str(p.relative_to(OUT)):build.sha256(p) for p in files},
  'frozen_v11_dependencies_verified':baseline['dependencies'],'frozen_v11_snapshots_verified':frozen['generation_snapshot_sha256'],
  'runtime':{'python':sys.version,**{n:importlib.metadata.version(n) for n in ['numpy','scipy','pandas','requests','threadpoolctl','anndata','h5py']}},
  'additional_dependency_sha256':{'record_v10_result.py':build.sha256(ROOT/'record_v10_result.py')},
  'v11_model_sha256':build.sha256(ROOT/'experiments/v11_consensus_coverage/analysis/reference_model.json'),
  'scope':'Research records only. No fitted replacement model, generated submission, upload or V11 numerical mutation.'})
 print('Research report and provenance written; V11 frozen inputs unchanged.',flush=True)
if __name__=='__main__':main()
