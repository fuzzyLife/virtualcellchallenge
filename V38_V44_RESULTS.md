# V38–V44: evidence and remaining gaps

Best published reference remains V24, overall **0.14299785521745675**. These raw public metrics use repeatedly examined biology and are not official anchored scores.

| Version | New question | Result |
|---|---|---|
| V38 | Response activity/direction/magnitude; K562 supervision | No development candidate passes; V24 retained |
| V39 | Independent KOLF experimental response features | No development candidate passes; V24 retained |
| V40 | Contrastive intervention-identity loss and shuffled-label control | No development candidate passes; V24 retained |
| V41 | Matched-model contrast and source-validation abstention | No development candidate passes; V24 retained |
| V42 | Smaller5%/10% identity calibration | No development candidate passes; V24 retained |
| V43 | Independent LINCS RNAi functional fingerprints | Fails matched shuffled-fingerprint mean control |
| V44 | V24 means, increased empirical residual variation | No development candidate passes; V24 retained |

## What the failures taught us

V38 gains were largely reproduced by removing the learned head. KOLF agreement is strongest for H1 and weak across other families. V40 improved source objectives but reduced generated-cell discrimination and coverage. V41 improved coverage while worsening other metrics; smaller V42 gains did not fix the tradeoff. Full empirical population transport repeatedly reduced false control-null DE calls but collapsed treatment fidelity. All failed outcomes and amendments are retained.

## New experimental fingerprints

Downloaded LINCS GSE106127 CGS PRIME RNAi consensus profiles, verified publisher SHA512s, and checked units. They are moderated z scores, not LFC or raw counts. RNAi seed/off-target effects remain. Quality filtering produces3,722 target fingerprints, including108/300 panel genes. PCA64 uses956 landmarks after removing all300 target outputs; each other intervention own landmark is zeroed. The features contain no held CRISPR responses. Source-specific bilinear ridge trains on2,360 K562 and559 RPE1 targets, plus eligible sources, excluding the complete held study from prediction. Missing fingerprints/outputs retain V24. A matched fingerprint permutation tests target-specific information.

| V43 arm | LFC squared error / V24 | Absolute error / V24 | Ideal pooled PDS |
|---|---:|---:|---:|
| fingerprint125 | 0.996487 | 0.998320 | 0.743687 |
| fingerprint250 | 0.995876 | 0.998639 | 0.730610 |
| shuffled250 | 0.995420 | 0.998309 | 0.726732 |

The genuine smaller blend improves ideal discrimination, but the shuffled control matches or improves mean errors. This bounded pilot does not show reliable useful RNAi transfer, so it is not promoted.

## V44: keep responses and change population dispersion

V24 combines25% expectation,65% Poisson counts,10% empirical residual before dual-moment fitting. V44 tests20% or35% empirical residual, compensating the expectation component; an adaptive arm uses 0.1+0.25*0.15/(0.15+abs(predicted LFC)). Strong effects retain more compact variation. The65% Poisson draw, exact libraries, control sampling and V24 response means are unchanged. Own-target emission keeps10%. Unlike prior full empirical transport, these arms retain compact sampling.

Five split-control nulls and20 treatment populations compare all six raw metrics, before locking a passing arm for24 additional count cases at two seeds. D3 remains unopened. Biological means were already used publicly, so count confirmation is not described as blind lineage validation.

| Development arm | PDS | FID | REACH | JAC | NMAE ↓ | Expression ratio ↓ |
|---|---:|---:|---:|---:|---:|---:|
| v24 | 0.791667 | 0.581703 | 0.194389 | 0.070027 | 0.912601 | 37.471936 |
| residual20 | 0.791667 | 0.573955 | 0.227109 | 0.068595 | 0.916782 | 37.448580 |
| residual35 | 0.805556 | 0.548126 | 0.232752 | 0.062733 | 0.912533 | 37.398750 |
| adaptive | 0.805556 | 0.545597 | 0.230151 | 0.063231 | 0.914251 | 37.399412 |

## Reproduction

Scripts are in the root. Use Python3.12 with requirements-v35-science.txt for scientific work, and .competition_env for official packaging. SCIENCE_PYTHON below denotes that pinned interpreter. Earlier V7/V30/V34/V37 documented source/cache workflows are prerequisites; source reuse requires matching hashes. Preserve existing receipts rather than overwriting experiments.

~~~bat
SCIENCE_PYTHON collect_v43_lincs.py
SCIENCE_PYTHON unpack_v43_lincs.py
SCIENCE_PYTHON prepare_v43_lincs.py
SCIENCE_PYTHON prepare_v43_fingerprints.py
SCIENCE_PYTHON test_v43_fingerprints.py
SCIENCE_PYTHON run_v43_fingerprints.py freeze
SCIENCE_PYTHON run_v43_fingerprints.py train
SCIENCE_PYTHON run_v43_fingerprints.py prepare
SCIENCE_PYTHON run_v43_fingerprints.py means

SCIENCE_PYTHON test_v44_response_dispersion.py
SCIENCE_PYTHON run_v44_response_dispersion.py freeze
SCIENCE_PYTHON run_v44_response_dispersion.py prepare
SCIENCE_PYTHON run_v44_response_dispersion.py means
SCIENCE_PYTHON run_v44_response_dispersion.py null
SCIENCE_PYTHON run_v44_parallel_counts.py development
SCIENCE_PYTHON run_v44_response_dispersion.py select
REM Continue only after a new arm passes development:
SCIENCE_PYTHON assess_v44_release.py lock
SCIENCE_PYTHON run_v44_parallel_counts.py confirmation
SCIENCE_PYTHON assess_v44_release.py assess
SCIENCE_PYTHON audit_v38_v44_results.py --version 44
SCIENCE_PYTHON audit_v44_raw_scores.py
REM Generate only after confirmation passes:
SCIENCE_PYTHON build_v44_submission.py fit
SCIENCE_PYTHON build_v44_submission.py preflight
SCIENCE_PYTHON build_v44_submission.py generate --workers 4
.competition_env/Scripts/python.exe package_v44_submission.py
SCIENCE_PYTHON audit_v44_delivery.py
~~~

Parallel context execution preserves numerical code/seeds/splits and isolates aggregate outputs, merging after all workers succeed. The execution amendment records that optimization after initial count outcomes. An initial object-string serialization failed before model fitting; its receipt is preserved and corrected Unicode arrays load without pickle. Large models/data/packages are ignored by Git.

## Primary sources

- [LINCS L1000](https://doi.org/10.1016/j.cell.2017.10.049)
- [RNAi versus CRISPR consensus/seed effects](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2003213)
- [Official GEO guide](https://clue.io/GEO-guide)
- [Publisher files/checksums](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106127/suppl/)
- [Cross-cell-line biological-response features](https://www.biorxiv.org/content/10.64898/2026.09.28.755146v1.full), a recent chemical-perturbation preprint; not evidence of genetic gains here.

Official wrapper validation only checks a nonempty member. Delivery additionally decompresses and compares every CSR value/index/pointer and all axes with the staged source. The final independent audit checks360,000 row libraries and900 groups. No automatic leaderboard upload or Git push.
