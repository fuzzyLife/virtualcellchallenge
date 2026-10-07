# Cross-modal foundation model review for the VCC project

Reviewed 2026-10-07. The useful idea is to improve biological representations using additional modalities and contexts, then learn a perturbation-response head on measured interventions. This is worth a controlled pilot. The supplied project does not currently provide a usable proteomics-pretrained checkpoint, and its protein reconstruction benchmark does not establish a gain on our CRISPRi task. V4 remains our production comparator; no submission was changed by this review.

## What was inspected

Read all 11 pages of the supplied `cross_modal_foundation_model/paper.pdf`, including supplementary methods, and rendered Figure 3 on page 6 to verify the bar labels. Inspected the reproduction's benchmark, normalization experiments, LFQ/TMT parsers, Phase 1 evaluator, educational transformer, dependencies, and local history. The audit hashes all 22 Python files and records function locations. Large benchmark arrays are inventoried by size; they were not loaded in full or treated as training labels for VCC.

The local Git remote is `https://github.com/animesh/cross_modal_foundation_model`, at `b4541490023df746b6d94375a9eaa7bd7bb1141b`. The history is:

- `f7718a5` (August 14): original author commit containing only the license.
- `0830655` (August 26): local protein baseline verification and PRIDE parsers.
- `eb473c3` (September 19): metadata reconstruction, normalization experiments, rank metrics and benchmark changes.
- `b454149` (September 30): exclude nine auxiliary/pilot LFQ samples absent from the formal metadata.

This is an extended reproduction, not the author's foundation-model training implementation. Some local files are untracked within the clone; `code_inventory.csv` identifies those explicitly. Their hashes, rather than the Git commit alone, pin the reviewed versions. The tracked clone was clean. The live folder contains more than 3,000 files and roughly 440 GiB of data; inventory totals are a point-in-time observation, not a frozen snapshot.

## What the paper actually does

Burq et al., *Single-cell foundation models benefit from cross-modal training: adding proteomics data beats parameter scaling*, bioRxiv DOI [10.64898/2026.08.14.744845](https://doi.org/10.64898/2026.08.14.744845), version posted August 19, 2026; supplied PDF SHA-256 is recorded in `review.json`.

The authors start with the RNA-trained Tahoe-x1 70M checkpoint, map proteins to gene tokens, quantile-bin abundances within each sample, and continue pretraining for one epoch on 48,843 proteomic samples from 440 studies. Raw and batch-corrected views provide augmentation. This is unpaired continued pretraining in a shared gene-token space, not a trained RNA-to-protein translator or a cell-mixture deconvolution method. Its outputs are gene/cell representations, not our required raw-count CRISPRi submission.

The paper reports gains over base Tahoe checkpoints on gene essentiality, pathway membership and tissue separation, and includes bulk-RNA continued-pretraining arms and multiple seeds. Those controls support investigating complementary data; they do not isolate every possible effect of compute, preprocessing or downstream evaluation. In particular, one epoch over the RNA and protein arms means different numbers of samples/updates.

Figure 3 uses leave-one-drug/context-pair-out cosine kNN to reconstruct a proteome from other measured perturbed pairs, then subtracts a matched control. Top-100 proteins are selected by the measured response for scoring. Approximate primary correlations, checked against the figure and paper text, are:

| Representation | Cardiomyocyte protein benchmark | Cancer-line protein benchmark |
|---|---:|---:|
| Base Tahoe-x1 70M | 0.542 | 0.746 |
| Proteomics continued pretraining | about 0.611 | 0.813 |
| HVG | about 0.665 | about 0.902 |
| PCA | about 0.674 | about 0.906 |

The negative annotations above the small gray bars are changes relative to base 70M, not negative HVG/PCA correlations. Some reproduction docstrings still misidentify these bars despite corrections in its README. Cross-modal pretraining improves its base model, but HVG/PCA are stronger on this particular benchmark. Paper Pearson correlations cannot be converted into expected VCC leaderboard scores.

## Why that evaluation cannot be copied into our submission

Our A/B/C matrices contain non-targeting controls only. We must predict 18,533 RNA counts for 400 cells per target/context across 300 CRISPRi interventions. The hidden treated profile is unavailable.

In `phase0/benchmark.py:345`, HVG queries `_median_centre(pair.r)`; PCA transforms `pair.r` at line 381. Here `pair.r` is the observed held-out treated proteome. The paper also describes embedding all available proteomic samples before predicting a held-out pair (pages 4-5). This makes the exercise reconstruction from an observed treated state. It is a different information setting from predicting an unseen intervention from controls and its identity. The authors' unpublished evaluation code could clarify details, but the local implementation definitely uses the treated query.

The neighbor pool at line 462 excludes only the current drug/context pair. It can retain other contexts for that same drug and other pairs marked held out from pretraining. Holding a drug out of pretraining does not hold every response to that drug out of the supervised neighbor pool. A leave-one-pair-out reconstruction score is therefore neither a fully unseen-drug predictor nor an unseen-context CRISPRi validation score.

VCC also evaluates direction, DE significance, effect magnitude, discrimination and expression error using generated cell distributions. A model can improve correlation or reconstruction error and still worsen DE fidelity if it creates weak wrong-sign effects or suppresses cell variance. We observed that with prior neural smoothing and the V5/V6 tradeoffs. V4's overall score (-0.050674) remains better than V6 (-0.056790).

## Specific code issues verified

These concern the local reproduction, not an assertion that the unavailable author implementation has identical defects. The new audit executes only named pure helpers on small synthetic arrays; it never invokes download/training entry points.

| Location | Finding | Consequence |
|---|---|---|
| `phase0/benchmark.py:177`, HVG at 339-346 | `_median_centre` replaces missing proteins with zero before the supposedly pairwise-complete similarity | A synthetic `[1, NaN, 3]` becomes `[-1, 0, 1]`; missing proteins now participate in similarity. |
| `phase0/benchmark.py:358-383` | PCA fits centered training profiles but transforms an uncentered held-out profile | Training and query use different preprocessing. |
| `phase0/benchmark.py:236`, pairwise kNN below it | Weighted `nansum` does not renormalize weights per detected protein | With equally weighted values `[NaN, 4]`, prediction is 2 rather than the valid-value weighted mean 4. |
| `phase0/benchmark.py:214` | Fewer than two valid features returns two values, while the caller unpacks four metrics | Sparse edge cases can raise an unpacking error. |
| `phase0/test_normalisation.py:128-161` | Batch means and HVG selection use all pairs, including evaluation pairs | This experiment has transductive preprocessing; it should not select our VCC model settings. |
| `phase1/run_phase1.py:282-283` | Hallmark logistic regression trains and scores on the same genes | Reported AUPRC is training performance, not held-out generalization. |
| `phase1/run_phase1.py:246-247` | Essentiality uses other cell lines' aggregate labels and evaluates the same gene identities | It tests transfer of shared gene essentiality labels; it does not evaluate a CRISPRi RNA response vector or unseen-gene generalization. |

The saved LFQ CSV has 59 scored pairs, versus 58 in the paper. Recomputed from that stored CSV, mean top-100 Pearson is HVG 0.556852, PCA 0.373533, perturbation mean 0.469130 and additive linear 0.462293. These are cached reproduction results, not new runs or proof that the September 30 parser produces them. Parser changes, pair identities and data processing must be reconciled before claiming replication. The ProTargetMiner manual channel map remains empty; no usable TMT benchmark result was established here.

Mass-spectrometry non-detections are encoded as NaN in this parser. Our RNA count zero is an observed count and must not automatically be discarded or converted to missingness. Preserve separate missing-value masks and the existing pairwise-valid rank strategy; avoid substituting protein median-centering for RNA count normalization.

## What is available to use now

Checked the public release manifests on October 7 and saved responses in `public_release_check.json`:

- [Author GitHub repository](https://github.com/tesorai/cross_modal_foundation_model): only `LICENSE` in the public root listing.
- [Trove1-70m release](https://huggingface.co/Tesorai/Trove1-70m/tree/main): only `.gitattributes` and `README.md` at commit `f1629213607e9e6a347138de35ee754a578e0738`; no listed model weights, tokenizer or configuration. Hugging Face's public API succeeded although the page failed through the browser search tool.
- No `.pt`, `.pth`, `.ckpt`, `.safetensors` or `.bin` checkpoint candidates were found in the local clone. The educational `minimal_llm_numpy_explainer.py` is a character-level NumPy text transformer, not a biological checkpoint or cross-modal pretraining pipeline.

There are useful precomputed **RNA-only Tahoe** gene embeddings in the local benchmark data. These are not Trove/proteomics-trained embeddings. Exact, unique symbol alignment using `obs.gene`, rather than its numeric row index, gives:

| Local mean embedding | Dimensions | Challenge genes covered | Challenge targets covered |
|---|---:|---:|---:|
| Tahoe-x1 70M | 512 | 13,980 / 18,533 | 235 / 300 |
| Tahoe-x1 1B | 2,048 | 13,980 / 18,533 | 235 / 300 |
| Tahoe-x1 3B | 2,560 | 13,980 / 18,533 | 235 / 300 |
| PCA | 15 | 13,982 / 18,533 | 236 / 300 |

Each file is hashed; all 300 target mappings are recorded. These cached gene representations are usable as optional features after verifying their upstream processing and reuse terms. They do not provide an encoder for new A/B/C cells or measured perturbation-response labels. Missing targets and output genes require explicit fallback; never silently fill missing embeddings and assume equivalent coverage.

## Recommended pilot, before another submission

1. **Use V4 as the fixed baseline.** Preserve the existing empirical control-cell generator. Compare all candidates on identical cell samples and random seeds, and report absolute raw metrics as well as aggregate score proxies. V6 did not beat V4, so it should not become the default teacher.
2. **Train on measured intervention contrasts.** Start with K562, CD4 and H1 response data already in this project. Align gene symbols, create matched perturbed-minus-control contrasts, and keep donor/batch information. Training targets are downstream RNA response vectors, not the observational 300-to-other-genes reconstruction task. A CRISPRi input-gene reduction fed through our old control decoder is not a learned causal response.
3. **Compare representations at small scale.** Use pairwise-valid Spearman/HVG/PCA baselines and the available 512-dimensional Tahoe gene vectors as a target prior. Start with strongly regularized ridge or a small low-rank residual head. Learn how target features and a control-derived context representation predict response deltas. Include a baseline without embeddings, and a rank/PCA context baseline with the same response head, so any gain can be attributed to the extra features. An engineered pooling of gene vectors is a new context representation to test; it is not the paper's released cell encoder.
4. **Validate both relevant transfer settings.** Hold out an entire source context for testing cross-context transfer, using only its controls at prediction time; fit normalization, feature selection, neighbors and response priors inside training folds. Separately hold out targets to measure new-target transfer. Known target responses in other source contexts may be used in the first setting, but must be clearly distinguished from unseen-target evaluation. Group related targets and source batches where possible. Do not use held-out treated expression to construct queries. H1's current batch split is useful for noise calibration but is insufficient by itself to demonstrate a new-context gain; the small number of public contexts limits precision.
5. **Gate adoption on VCC-relevant behavior.** Evaluate generated counts with the existing rank-DE implementation and the panel gene exclusions. Require fidelity, DE overlap and expression error to improve consistently over V4 on held-out data, while retaining discrimination, reach and cell dispersion. Report target-level variation and bootstrap uncertainty, not only mean correlations. The nine larger H1 target groups are a small development sample and have already been inspected; success there is not independent confirmation of hidden A/B/C accuracy.
6. **Test genuine cross-modal weights when released.** Freeze the encoder first and train only the response head. Compare against RNA-only Tahoe and classical representations with the same data, splits, head and generator. A full 440-study proteomics curation/pretraining effort is a later option if this inexpensive pilot supports the representation hypothesis and the original checkpoint remains unavailable.

The most promising near-term use is an embedding-assisted, context-aware response residual trained on real perturbations, with conservative fallback to V4. A richer prior might improve response signs or target coverage, but neither improvement has yet been measured. Directly copying the paper's kNN reconstruction or making a larger control-expression neural network would not address the missing interventional evidence.

## Reproduce this review

Run in Windows Command Prompt from the project root:

```bat
cd /d L:\promec\Animesh\Download\Virtual
.\.competition_env\Scripts\python.exe review_cross_modal_model.py
```

Optional current public-release verification:

```bat
.\.competition_env\Scripts\python.exe review_cross_modal_model.py --check-online
```

Outputs are in `analysis/cross_modal_review/`: code hashes/function locations, input provenance, embedding coverage/target maps, stored-baseline summaries, synthetic code probes and the dated public release check. The extracted full paper text and rendered intermediate figure are ignored. Dependencies are the existing project environment plus `pypdf`; no new checkpoint, training corpus, pretrained model inference or leaderboard upload is required. The synthetic checks verify code behavior, not improved competition performance.
