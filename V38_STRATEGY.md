# V38: response activity and conditional magnitude

## What is new

V24 remains the published reference (overall 0.1429978552). V36 and the seven
V37 ensembles did not establish an improvement. V29 preserved empirical cell
variation and passed control-null tests, but lost treatment direction yield.
Therefore this experiment changes the **supervised response objective**, and
only then tests whether empirical variation transport becomes useful.

Three fixed heads/corrections:
1. A natural-prevalence three-class classifier (down / inactive / up) and
   separate positive/negative conditional magnitude regressors. The posterior
   expected signed magnitude corrects the pooled response by one quarter.
2. A residual regressor with eightfold weight on active responses, again at
   quarter strength. Unlike ordinary mean-square fits, this objective does not
   let inactive genes dominate the loss.
3. A half-strength selective correction only where activity probability is at
   least 0.20, directional confidence at least 0.8, and the learned effect
   agrees with and strengthens the pooled effect.

All corrections are bounded to one log2 unit before their fixed coefficient;
the own-target response is zero in the transfer head. The count generator
applies the separate established 90% own-target assumption once.

## Evidence informing the hypothesis

[scDEBART](https://proceedings.mlr.press/v306/sung26a.html) learns expression
changes conditioned on basal expression rather than reconstructing static
co-expression. Its authors report benefits from response-specific pretraining.
V38 is an independent small supervised experiment, not an implementation of
scDEBART or a use of its pretrained weights.

[The official scoring definition](https://raw.githubusercontent.com/ArcInstitute/cell-eval2/main/docs/vcc2026_metrics/vcc2026-metrics.md)
uses Wilcoxon-derived response statistics as well as pooled-count profiles.
Matching expression means does not by itself match significance or direction.
The scorer document and the upload CLI have different scope: packaging follows
the competition's 400 cells per target/context and excludes added controls.
Private anchor scaling is not reproduced by our public diagnostics.

## Data, exclusions and model

The checksum-verified V7 cache contains 1,462 recipient/target rows, 300 target
genes and 18,533 output genes. K562 is now also a supervised recipient; the
recent V35 table used its experiments as retrieval donors only.

Recipient training uses at most 256 uniformly sampled valid output genes per
row; there is no enrichment of positive labels or case-control probability
distortion. Equal biological-study weights prevent the numerous CD4/Orion
rows from drowning out H1. Orion groups HCT116 and HEK293T into one source
study and includes their disagreement in uncertainty.

For every outer fold, exclude **all** H1, Orion, or CD4 outcomes from fitting.
For each fitting recipient, also exclude its entire study from experimental
donor features. Features contain basal gene/target expression and ranks,
STRING coordinates, supported donor effect/uncertainty, agreement and control
similarity. Missing donors are NaN features, never biological zero votes.
An explicit mutation test changes excluded H1/K562 outcomes and verifies
unchanged donor features. All occurrences of hash-zero targets are excluded
from fitting heads and used for source validation; allowed other-study
experimental responses to these interventions remain retrieval features.

Activity labels require absolute log2 effect >=0.25 and at least twice its
cached standard error. The cache combines source pipelines; these are
**diagnostic effect/uncertainty labels, not competition Wilcoxon DE labels**.
Generated cells are evaluated with actual Wilcoxon/BH tests.

Each head uses 120 histogram boosting iterations, 15 leaves, minimum leaf
size 80, L2=10, learning rate 0.06, fixed seed 38038. No held-query outcomes
enter inference. Model files, inputs and scientific dependencies are hashed.

## Fixed tests and selection

The frozen executable plan is
`experiments/v38_response_hurdle/analysis/plan.json`.

All three models first undergo the same 44 previously spent public mean
comparisons. The lowest equal-study mean-square-error ratio among candidates
with MSE/MAE within 1% of V24 proceeds to count tests. If none meets that bound,
the best candidate is diagnostic only and cannot qualify automatically.

Compare V24, that selected mean head with the established generator, and
the same head with the V29 empirical anchored transport. Test the original
first four cases per context, all six public raw diagnostics, and five
three-way control-null splits. Templates belong to fitting controls; actual,
reference and fitting populations are mutually disjoint. Additional 24 count
cases and a second generation seed are reserved for confirmation of a chosen
candidate. Public mean labels/source biology have already been used: this
is not a new blind lineage test.

Frozen count guards: PDS within 0.01, fidelity within 0.01, reach within 0.005,
Jaccard within 0.001, NMAE within 1%, expression error within 2%; all studies
resolved; null DE calls no worse in every context; positive average relative
gain over the six raw metrics. These are development decision rules, not a
prediction of private score or statistical significance.

## Reproduce

Use the scientific Python environment with the pinned V35 science dependencies.
On this machine:

```cmd
cd /d L:\promec\Animesh\Download\Virtual
set PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe
%PY% test_v38_response_hurdle.py
%PY% run_v38_response_hurdle.py freeze
%PY% run_v38_response_hurdle.py train
%PY% run_v38_response_hurdle.py prepare
%PY% run_v38_response_hurdle.py means
%PY% run_v38_response_hurdle.py null
%PY% run_v38_response_hurdle.py development
%PY% run_v38_response_hurdle.py select
```

Frozen outputs are write-protected by existence/hash checks. For a fresh replay,
copy the scripts and set V38's output namespace to a separate replay folder;
do not overwrite the frozen plan, models or receipts. Neither this runner nor
packaging uploads a submission.


## Binning compatibility amendment

The first fit stopped before producing a model: the installed sklearn/NumPy
combination could not bin an entirely missing donor column. The frozen initial
plan and code are retained. The adapter removes training-constant/all-missing
columns separately per head; no validation rows select columns. A small
constant/missing-column regression check and an exact joblib round-trip passed.
The first adapter launch was stopped before fitting to use an importable
class name for model serialization. Use this training command instead:

```cmd
%PY% -c "import v38_binning_adapter as a; a.main()"
```

The amendment and adapter hash are recorded in analysis/binning_amendment.json.
