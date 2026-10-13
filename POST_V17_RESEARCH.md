# Research after V17: methods we have genuinely not tested

Research and code/asset audit: **2026-10-13**. Best published project artifact remains V16, overall 0.141613857. This work downloads small author code, metadata and a gene vocabulary; it does not load model weights, train models, generate predictions or alter any submission.

## Main recommendation

The strongest next experiment is to introduce an independently trained **population-effect expert with an explicit DE-confidence head**, then test a **pretrained population generator with paired synthetic controls**. Specifically: evaluate Arc PIE first, and Stack second. Both have public checkpoints. Their architectures and inputs differ substantially from our source-pooling variants, control-reconstruction networks and small response heads. Whether their prediction errors are complementary remains untested.

Our present model often distinguishes perturbations without recovering their exact DE sets. V15 experts had residual-error correlations around 0.99–1, so voting among them added little. V17 changed a single response-strength variance while retaining the V16 average effects, supported only 85 targets, and failed its shuffled-target guard. The next expert should bring new information and predict which genes respond, not merely change the spread around the same response vector.

## 1. PIE: a new response expert, available now

Arc announced [PIE on October 5, 2026](https://arcinstitute.org/news/pie); the [preprint](https://doi.org/10.64898/2026.10.02.756297), [author code](https://github.com/ArcInstitute/pie) and [model collection](https://huggingface.co/collections/arcinstitute/pie) are public. PIE combines biological knowledge, control expression and training-response evidence to predict population-level effects. Its [checkpoint card](https://huggingface.co/arcinstitute/PIE_replogle_xdataset) specifies three outputs: `p_de`, `lfc_pred` and `delta_p_pred`. These are not 400 generated cells.

**What is new for us:** joint DE classification and effect prediction using pretrained biological knowledge. We have not evaluated a PIE checkpoint. Our Reactome imputation and Tahoe response heads do not establish whether these richer features and learned heads are useful.

Our exact metadata audit finds:

| Published asset | Competition output genes present | Panel targets present |
|---|---:|---:|
| PIE gene-description embeddings | 18,262 / 18,533 | 300 / 300 |
| PIE perturbation-specific text | 10,710 / 18,533 | 155 / 300 |
| PIE ESM2 protein embeddings | 17,876 / 18,533 | 299 / 300 |

These are feature-availability counts, not proof that a checkpoint predicts every coordinate reliably. Missing descriptions must remain missing; gene-description coverage does not imply complete perturbation-text coverage.

**Critical evaluation issue:** the cross-dataset checkpoint trained on Tahoe, Jiang, ARC VCC25 and Orion. Our H1/HCT116/HEK293T count tests use ARC VCC25/Orion biology and may include its training cells. Those tests cannot establish unseen-family performance for that checkpoint. Start with the [K562-held-out within-dataset checkpoint](https://huggingface.co/arcinstitute/PIE_replogle_wdataset), trained on HepG2/Jurkat/RPE1 instead. We checked its exact dataset revision `20c9faef76fc96fdc809871340bd93499413dfe7`: none of its 1,857 perturbation names is in our 300-target panel. All our targets are therefore unseen perturbation identifiers for this checkpoint, not guaranteed predictable ones.

The within-dataset training output vocabulary has 6,642 genes, of which 6,290 overlap ours. First test this covered vocabulary while retaining V16 elsewhere. Extension to additional genes through text queries is a separate extrapolation experiment. Do not claim 18,262 trained output genes just because descriptions exist.

`p_de` is a model probability, not a Wilcoxon/BH q-value. It must be validated at our 400-cell setting. Different original source sample sizes and DE-label construction affect its calibration. Unknown A/B/C identities should use control-expression inputs and missing context metadata, rather than guessed Cellosaurus identities. Verify that source masking works in actual inference before relying on it.

**Feasibility:** one within-dataset checkpoint is 949,595,187 bytes; the cross-dataset checkpoint is 966,412,013 bytes. The code requires Python 3.12, while our current environment is Python 3.14.3. Full configured assets are much larger: the published NCBI-text embedding array alone is 43,182,215,296 bytes. Inspect/download a bounded set of sources or test source masking before downloading the whole stack. Original training-response assets remain necessary for the official inference path; the checkpoint alone is insufficient.

## 2. Stack: transfer a population effect using cells as the prompt

[Stack](https://arcinstitute.org/tools/stack) uses cell populations as examples at inference. The [author generation tutorial](https://github.com/ArcInstitute/stack/blob/main/notebooks/tutorial-predict.ipynb) transfers drug effects from T cells to control myeloid/B cells. Its [public model card](https://huggingface.co/arcinstitute/Stack-Large) describes a 217M-parameter model. This is an attractive match to anonymous A/B/C contexts because it can use actual control cells without requiring their cell-line names. Genetic-knockdown transfer still needs testing; the drug tutorial is not evidence of success on our CRISPRi panel.

**What is new for us:** inference with a pretrained population model and observed source treatment cells as prompts. We have not run Stack. V10's integer moment transport, V14's two-feature context interaction and V17's scalar variance do not learn this conditional population map.

The tutorial explicitly uses a generated-control reference to remove DE artifacts. That matters because V17 recorded 1,946.5 generated-NTC DE calls against independent real controls, versus zero calls for the actual diagnostic controls. This supports investigating generator baseline mismatch; it does not prove that all discrepancies arise from one cause.

For each source and target, generate both:

- source NTC prompt + query-context controls -> generated query baseline;
- source treated prompt + the same query controls -> generated query treatment.

Estimate a treatment residual relative to that generated baseline. Then calibrate the residual onto the **real** query baseline, because the competition scorer uses its real controls. Scoring only generated treatment versus generated control would conceal a submission-time mismatch. The residual-to-count conversion is an experiment in its own right; paired seeds do not automatically cancel every nonlinear count bias.

Our non-executing vocabulary audit finds 15,012 model genes, including 13,078 competition genes and 213 panel targets. The model must therefore retain a validated fallback for 5,455 output genes. A missing target coordinate need not make the downstream prompt unusable, because Stack conditions on the observed population, but that is a testable possibility rather than established coverage.

The aligned checkpoint is 2,613,863,242 bytes. Our RTX 2070 SUPER has 8 GiB and working CUDA, but no native bfloat16. A batch-one memory/precision smoke test in an isolated compatible environment is required before promising a full run. The upstream alignment helper zero-pads missing genes; those placeholders must not become submission predictions. The inspected prediction path returns sampled negative-binomial counts, and its gene alignment/library handling must be audited before combining them with fallback counts.

Stack's large public pretraining corpus may overlap public test datasets. Its cell-ID/dataset overlap is not audited here. Label retrospective results accordingly and seek an independently audited source for strong generalization claims.

## 3. Covariance-response transfer: a small, mechanistically different CPU experiment

The [CIPHER preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC12363937/) proposes a linear-response relation between control covariance and perturbation-induced mean shifts. Its [current forward code](https://github.com/GoyalLab/CIPHER/blob/main/cipher/forward.py) fits a response amplitude using observed treatment shifts on training genes; optional gene holdouts assess the remaining genes. That is not the competition's no-treated-query-data setting. The full older preprint and current implementation were inspected; the September 2026 revision's abstract was accessible, while direct full-text access was unavailable.

A deployable adaptation would instead estimate a target's force vector from **other source contexts**, then apply the query control covariance. In a training-only shared program space:

```text
source_force[target] = regularized_inverse(source_control_covariance)
                       @ source_observed_mean_shift[target]
query_mean_shift[target] = query_control_covariance @ pooled_source_force[target]
```

This can produce context-specific directions even when control means are similar. It differs from V14's regressions on two control-mean rank components and V17's response-amplitude variance. It also differs from treating our 300-to-other-gene reconstruction network as a causal simulator.

Use 8–32 shared programs, depth/batch-adjusted controls, masked measured genes, shrinkage of unstable covariance directions, and whole-family exclusion. Do not independently whiten every context to identity and thereby erase the covariance differences being tested. Compare real covariance with diagonal covariance, shuffled program labels and ordinary conserved response transfer. No A/B/C treated labels may estimate the force.

A related population model exponentially reweights real query cells using learned source forces. This could change relative state occupancy rather than add generic noise. [NubOT](https://github.com/FrederikeLuebeck/nubot) provides a precedent for combining transport and reweighting; it does not establish that our snapshots identify birth/death rates. First test whether held-source treated distributions are even representable by reweighting control states. If treatments leave control support, reweighting alone cannot solve the task.

## 4. Conditional flow or copula models: useful second-stage options

[CellFlow's author documentation](https://cellflow.readthedocs.io/en/latest/notebooks/100_pbmc.html) distinguishes deterministic OT flow matching from stochastic `GENOT`, with sample/context representations. [Conditional Monge Gap](https://www.nature.com/articles/s42256-026-01242-8) learns one conditional transport model rather than one map per condition. These are genuinely untrained model classes here, but their demonstrated drug/cytokine settings do not ensure unseen-lineage CRISPRi performance. Both need suitable cell-level source data and a tested count decoder.

[scDesign3](https://pmc.ncbi.nlm.nih.gov/articles/PMC11182337/) separates state-conditioned marginals from cross-gene dependence using copulas. We have tested unconditional gamma/Poisson counts, five-quantile transfer and a single latent response strength; we have not fitted a full state-conditioned marginal/copula model. Such a simulator can model variation, but cannot determine an unknown treatment response without a separately validated transfer model.

These approaches are best considered after testing the new response expert or demonstrating a remaining population-specific error budget. Training a larger generator before verifying that its response direction is correct repeats a known project failure mode.

## Data priorities from this review

1. Use existing verified H1/HCT116/HEK293T cells for the first clean PIE-checkpoint pilot. A/B/C controls suffice to build inference queries; their hidden treated data are unavailable and unnecessary for setup.
2. Investigate the CIPHER release's donor-resolved Marson/CD4 files. Its [pinned manifest](https://github.com/GoyalLab/CIPHER/blob/main/resources/zenodo_manifest.csv) lists D1 resting and stimulated files around 1.4–2.3 GB. These could add cell-level primary-T-cell validation to our currently summary-only CD4 evidence. Check actual raw-count layers, target overlap, donor/batch provenance, library axes and duplication before claiming new independent biology. They may repackage the same study we already use.
3. Selectively extract the KOLF atlas's QC-filtered cell-level matrix for matched targets/NTCs. V17 established 273 target names in its browser, not 273 powered raw-count groups. Prefer unselected QC-filtered cells over only strong-perturbation cells. The full current file is 189 GB, so begin with metadata and bounded reads.
4. Do not import the new manifest's fibroblast/RPE1 TF activation screens as if they were CRISPRi. [Southard et al.](https://pubmed.ncbi.nlm.nih.gov/40770575/) used CRISPRa; sign reversal is not a validated conversion to knockdown effects.

## What was completed and how to reproduce

- Inspected existing model implementations and failure reports; recommendations above distinguish implementation from proposal.
- Pinned six author repositories and four model repositories; archived selected source, model cards, licenses and asset metadata. No archived author code was executed.
- Verified 50 research inputs; safely interpreted vocabulary literals without executing pickle globals.
- Audited exact feature coverage, checkpoint sizes, exact PIE training-target overlap and current GPU/runtime compatibility. Reader/schema/coverage checks passed.
- Prepared the staged pilot in [POST_V17_NEXT_EXPERIMENT.md](POST_V17_NEXT_EXPERIMENT.md). No biological performance test of PIE, Stack, CellFlow, CIPHER or CMonge is claimed.

Root commands:

```bat
.\.competition_env\Scripts\python.exe collect_post_v17_research.py
.\.competition_env\Scripts\python.exe audit_post_v17_feasibility.py
.\.competition_env\Scripts\python.exe test_post_v17_research.py
```

Use `collect_post_v17_research.py --fetch-missing` to retrieve missing small immutable inputs; historical API snapshots come from Git. Exact URLs/revisions/hashes are in `analysis/post_v17_research/source_manifest.json`; measured coverage and runtime evidence are in `feasibility_audit.json` and `feature_coverage_summary.csv`. Model licenses and acceptable-use terms are archived where inspected; code and model-output licenses are distinct. Prior V16/V17 numerical files and artifacts remain unchanged.
