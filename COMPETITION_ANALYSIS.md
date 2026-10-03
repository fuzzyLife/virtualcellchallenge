# Virtual Cell Challenge 2026: project analysis

This project is the validation-round bundle for the 2026 Virtual Cell Challenge,
panel `vcc2026-val-1`. The root script `build_submission.py` analyses the controls,
transfers public genome-wide K562 CRISPR responses to the anonymous contexts, generates
raw single-cell count predictions, and runs the official validator and packager.

## Completed submission

**Upload `prediction.vcc`**, the official VCC package, rather than either reference
data or per-context control files. Its size is 3,306,905,600 bytes (3.31 GB).

- Official `vcc prep`: passed; target sets verified and counts preserved.
- Shape: 360,000 cells by 18,533 genes; 120,000 cells per context.
- Groups: 900, with exactly 400 cells in each context/target group.
- Sparse stored entries: 2,083,470,063, below the 4,750,000,000 cap.
- Largest predicted cell total: 52,457, below the 1,000,000-count cap.
- All stored values: positive whole raw counts; no explicit zeros.
- Gene order: exact match to the official CSV; no reordered genes.
- SHA-256: `d236c64869488827f54012a0e10de68aad05c6d97142674a561969483618371e`.
- Tests: 5 passed, including a 200-seed sampling regression check.
- Upload status: published, entry `QQnFEBGlOk2pFFwkgtJL` (user-provided CLI output).
- Rank at publication: 1140; overall score: -0.276289857599391.
- Full score record: `analysis/leaderboard_result.json`; interpretation and commands: `README.md`.

`prediction.h5ad` (2.47 GB) is the full intermediate prediction. The packaged file
has been checked with the official CLI's archive validator and its metadata matches
the independently scanned H5AD. Full evidence is in `analysis/submission_ready.json`
and `analysis/official_prep_stdout.json`.

## What the competition requires

| Requirement | This panel |
|---|---|
| Task | Zero-shot CRISPRi response prediction in unseen cellular contexts |
| Contexts | A, B, C; retain the labels on their original control files |
| Available inputs | 18,400 non-targeting control cells per context |
| Control guides | 46 guides, 400 cells per guide |
| Perturbations | Exactly the 300 symbols in `pert_counts.csv`, in every context |
| Cells per perturbation/context | Exactly 400 |
| Total predicted cells | 360,000; 900 context/target groups |
| Gene axis | All 18,533 symbols, in the exact order of `gene_names.csv` |
| Expression | Finite, nonnegative integer raw counts; no log normalization |
| Library-size cap | 1,000,000 counts per cell |
| Density cap | At most 4,750,000,000 stored entries; use sparse CSR |
| Control rows | Exclude all rows labeled `non-targeting` |
| Upload artifact | One `prediction.vcc` containing all three contexts |
| Packaging | Official `vcc prep`, with both reference CSV files |

The local gene list **has a `gene_name` header**. Treating it as headerless would
incorrectly add a feature named `gene_name`. `pert_counts.csv` contains only
`target_gene`; the 400-cell requirement comes from the manifest, not that CSV.

The six evaluation components measure perturbation discrimination, expression
accuracy, differential-expression effect-size accuracy, direction fidelity,
direction reach, and differential-expression overlap. Scores are scaled separately
for each context and averaged across the six metrics and three contexts. The
perturbed gene itself is excluded from every metric. Merely suppressing that gene
therefore does not solve the scored prediction task.

The validation leaderboard is distinct from the final competition result. The
provided documentation schedules the final inputs for October 22, 2026 and the
submission deadline for November 5, 2026. Final inputs use different cell contexts
D/E/F and a new perturbation panel; these validation predictions are not a final-round
submission. The script derives the context labels and counts from the manifest.

## Review of the project files

All 146 original files were inventoried with paths, sizes, SHA-256 hashes and their
role in `analysis/file_inventory.json`.

- `manifest.json`: authoritative season, partition, panel, shapes, labels and counts.
- `gene_names.csv`: ordered feature axis; 18,533 unique gene names.
- `pert_counts.csv`: 300 unique prediction targets, all present in the feature axis.
- `context_A.h5ad`, `context_B.h5ad`, `context_C.h5ad`: CSR raw count matrices;
  each has only non-targeting controls and the correct gene order. Every stored
  input count was checked for finiteness, nonnegativity and integrality.
- `downloads.html`: detailed bundle and submission specifications.
- `Virtual Cell Challenge.html`: registration, workflow and timeline information.
- `leaders.html`: saved validation leaderboard snapshot; not training labels.
- `FAQ.html`: saved FAQ topics; collapsed answers were not captured as visible text.
- The saved Cell article: explains the 2026 zero-shot, multi-context benchmark.
- `2025.09.12.675911v1.full.pdf`: a 61-page paper titled *Generative design of novel
  bacteriophages with genome language models*. It is unrelated to this prediction
  task and was excluded. Its text was extracted to establish its subject.
- The five `*_files` directories: saved browser JavaScript, styles, fonts/graphics
  and article assets. Their files were inventoried and text assets scanned; they
  do not contain a competition model, labeled response dataset, or model weights.

No challenge perturbation-response training set is supplied. All controls contain
`target_gene = non-targeting`, and their H5AD `layers` and `uns` are empty. The
withheld ground truth cannot be reconstructed from the manifest's cell counts.

## Control-data findings

| Context | Cells | Stored nonzeros | Median library counts | Median detected genes | Zero fraction |
|---|---:|---:|---:|---:|---:|
| A | 18,400 | 109,910,478 | 20,109 | 6,147 | 67.77% |
| B | 18,400 | 101,546,777 | 19,946 | 5,756 | 70.22% |
| C | 18,400 | 108,213,983 | 20,034 | 6,006 | 68.27% |

All 300 targets have nonzero expression somewhere in every context. Correlations
between basal log1p mean-CPM profiles are 0.7801 for A/B, 0.7478 for A/C, and
0.8409 for B/C. Contexts therefore have distinct basal profiles and must remain
separate throughout prediction.

The analysis includes per-gene raw mean, standard deviation, detection fraction,
and mean CPM for every context; target-only summaries; full library/detected-gene
quantiles; guide counts; and context-profile correlations. CPM is used only to
analyse and compare profiles. The uploaded prediction remains raw counts.

## Model and provenance

`download_reference.py` downloads the original authors' **K562 genome-wide raw
pseudobulk** release from Replogle et al. (2022), Figshare article 20029387, file
35774443 (`K562_gwps_raw_bulk_01.h5ad`). The 374,587,922-byte download was verified
against the publisher's MD5, `4570b53c9d62ff6df281e622f0350060`. The release is
CC-BY-4.0. Original file metadata and a provenance record are saved under `reference/`.
This file supplies pooled expression **means** from a public experiment, not hidden
challenge responses. The 65.8 GB raw single-cell file was not needed.

The reference has 11,258 profiles and 8,248 feature columns. Its gene symbols use
an older HDF5 categorical encoding in `var.gene_name`; its variable index contains
Ensembl IDs. The script decodes the symbols and sums duplicate symbols rather than
mistaking the IDs for the challenge gene axis. Construct names embed both the gene
and promoter/transcript design; these are parsed into gene symbols before pooling.

A first public essential-gene reference was checked and rejected because its target
panel was disjoint from the competition panel. It is not used for prediction.
`reference/rejected_essential_reference.json` records that coverage finding.

The retained genome-wide reference covers **260 of 300 targets** after requiring
at least 30 cells and a measured source knockdown of at least 50%. It covers
**7,681 challenge genes**. The other 40 targets and genes outside the source panel
receive a conservative control-derived downstream fallback. Their on-target counts
are still reduced according to the stated knockdown assumption.

The model:

1. Uses 503 source non-targeting core-control profiles, representing 75,094 cells.
   Pools effective source constructs for each gene with cell-count weights.
2. Aligns shared symbols and normalizes source means over that common gene axis.
   Computes log2 response ratios with a 1-CPM pseudocount.
3. Soft-thresholds changes using a Poisson sampling-noise approximation, clips
   extreme log2 effects at +/-3, and shrinks low-cell-count estimates. This is a
   noise-control assumption, **not** a formal differential-expression test.
4. Transfers 30% of each estimated log2 effect into each challenge context. The
   fixed 0.3 shrinkage is an explicit conservative modeling assumption; it is not
   fitted to challenge data or claimed to be validated across cell contexts.
   Each context retains its own basal expression and cell-distribution inputs.
5. Samples 400 distinct control-profile inputs per target, balanced over all
   46 guides (8 or 9 cells each) and stratified across each guide's library sizes.
   Integer quantile boundaries keep the sampling bins disjoint.
6. Simulates changed counts: binomial thinning for decreases, Poisson additions
   for increases. Uses 90% on-target knockdown as an assumption, not a fitted
   challenge measurement. All prediction rows receive new target/context labels.
7. Writes CSR blocks incrementally with 64-bit global pointers, removes zeros,
   and checks every count, library total, group, label and gene position.
8. Runs official `vcc-cli==0.2.2` validation and packaging. Only a successful
   package is renamed `prediction.vcc`; unfinished output has a `.partial` name.

`analysis/reference_cross_validation.csv` is an **exploratory construct-agreement
diagnostic**, despite that legacy filename. It compares different promoter or
transcript-targeting constructs for the same source gene. They are different
interventions, not independent technical replicates, so these values do not validate
cross-context transfer and are not used to tune the final transfer setting. The
weak agreement is an uncertainty about the source's transferable downstream signal.

## Validation and limitations

Five focused tests passed, including a 200-seed check against duplicate quantile-bin
samples and a synthetic H5AD-to-VCC roundtrip through the official packager. The
first full run exposed an integer-bin boundary overlap; this was fixed before the
complete prediction was rebuilt. `analysis/test_results.txt` records the checks.
Full prediction checks and the official result are recorded in
`analysis/prediction_validation.json`, `official_prep_stdout.json`, and
`submission_ready.json` when generation and packaging finish.

No local competition score is available because challenge responses are withheld.
There is no matched public replicate or challenge holdout with which to validate
performance in A/B/C. K562's tissue and sequencing assay differ from those unseen
contexts; transfer can fail. Coverage is partial, down/up effects are conservative,
and the simulator cannot induce a gene absent from all context controls. A Poisson
noise approximation does not model all biological or technical overdispersion.

The artifact is an experimentally informed statistical baseline. Local checks did not establish predictive performance. The first upload has now
been scored: rank 1140 at publication and overall -0.276289857599391, below the
competition baseline, principally because of direction fidelity (-1.6653311841128977).
No local hidden-response score or matched control-only submission is available.
These aggregate scores do not identify whether incorrect directions or insufficient
DE calls explain the fidelity loss. See `analysis/leaderboard_result.json`.

## Reproduce and upload

From this project root, using PowerShell:

```powershell
# The project-local environment has already been prepared.
& .\.competition_env\Scripts\python.exe .\test_submission.py
& .\.competition_env\Scripts\python.exe .\download_reference.py
& .\.competition_env\Scripts\python.exe .\build_submission.py

# Rebuild intentionally, keeping the same default random seed:
& .\.competition_env\Scripts\python.exe .\build_submission.py --force

# Repackage an existing full prediction without rerunning the model:
& .\.competition_env\Scripts\python.exe .\build_submission.py --package-only --force

# Authenticate interactively, then upload when you decide to submit:
& .\.competition_env\Scripts\python.exe -m vcc login
& .\.competition_env\Scripts\python.exe -m vcc submit .\prediction.vcc -m 'Public CRISPR transfer' --wait
```

`build_submission.py` packages locally; it does not upload or consume a daily
submission slot. API credentials are needed only for the final login/upload. Keep
keys out of scripts and version control. `.competition_env` is a local environment
and is not portable; on another machine install `requirements-competition.txt` into
a Python 3.14 environment. Allow tens of GB of temporary disk space and roughly
40-60 GB of memory for the official full-matrix packaging step.

## Sources

- Saved local competition and dataset pages, manifest, and reference CSV files.
- [Official CLI guide and submission requirements](https://vcc-cli-wiki.virtualcellchallenge.org/).
- [Official CLI command reference](https://vcc-cli-wiki.virtualcellchallenge.org/cli-reference).
- [Official evaluation page](https://virtualcellchallenge.org/evaluation).
- [Scoring implementation](https://github.com/ArcInstitute/cell-eval2).
- [Original-author public reference and CC-BY attribution](https://plus.figshare.com/articles/dataset/20029387).
