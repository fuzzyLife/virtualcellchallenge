# New response-learning tests: V38, V39 and V40

The task is still zero-shot context transfer. A new architecture is useful only
if it preserves perturbation identity, predicts credible effect sizes and
produces appropriate count populations. All experiments here reuse public
development biology; none proves a private leaderboard improvement.

## What we learned and tested

1. **V38: activity/direction and conditional magnitude.** Fit separate down /
   inactive / up classification and conditional effect-size heads, plus a
   response-weighted residual regression. Include K562 as an additional
   supervised recipient. Whole H1, Orion and CD4 studies are held out in turn;
   recipient studies are additionally excluded from their donor features.
   Activity uses effect/SE diagnostics, not heterogeneous cached q-values.
2. **V38 ablation: remove the learned head.** Mean-error gains largely came
   from attenuating the pooled response. A zero learned head performed slightly
   better. This prevents attributing the gain to successful activity learning.
3. **V39: add previously unused experimental evidence.** Download the pinned
   all-perturbation KOLF browser export, verify complete Git blob identities,
   and extract 273 targets x18,106 unambiguous output genes. Add seven auxiliary
   features to exactly the same V38 recipe. No changes to labels, folds, gene
   samples, coefficients or count generator.
4. **Empirical variation transport.** Test the V29 empirical anchor with the
   new mean heads. Corrected nulls improve, but held-treatment fidelity again
   collapses. Withdraw the emitter, preserving all partial evidence; finish
   only the unchanged standard-generator candidates.
5. **V40: train for response identity.** Fit small, positive per-gene gains for
   positive/negative baseline effects using a within-recipient contrastive
   objective, a mean-error auxiliary and regularization. Gains are bounded
   0.75..1.25; zeros stay zero and encoded effect signs cannot flip. Compare
   0.10/0.25 gain bounds. This changes the objective without inventing a new
   perturbation response from unperturbed correlations.

The separately fitted V40 negative control shuffles **contrastive positive
labels only**; its mean-error auxiliary remains correctly paired. It isolates
the value of the contrastive objective, rather than being a completely
outcome-shuffled model.

## Literature and applicability

[scDEBART](https://proceedings.mlr.press/v306/sung26a.html) motivates learning
expression changes conditioned on basal expression rather than reconstructing
static expression. Our V38 is an independent small supervised model, not its
architecture or checkpoint.

[PerturbCellRL](https://arxiv.org/html/2606.27752v1) motivates checking biological
response alignment beyond distribution matching. V40 directly trains a
differentiable response-identity objective; it is not an RL or flow-matching
reproduction.

The [KOLF atlas](https://www.nature.com/articles/s41587-026-03199-w) provides
independent genetic-intervention evidence in an iPSC context. Its
[author browser export code](https://github.com/y-doctor/KOLF2.1J_Perturbation_Cell_Atlas/blob/8f84531a784f24b14d7e2c9b80c375d7e78fd46c/psp/notebooks/mde_webapp/webapp_build/build_panel_all_perts.py)
writes int8 mean NTC z scores in log-expression space at scale0.04. These are
**not counts, log fold changes, DE probabilities or single-cell variation
examples**. We downloaded about392.5MB of browser binaries, not the189GB raw
AnnData or a strong-perturbation-only cell subset. Missing entries remain NaN;
quantization saturation has a separate flag. Every matching target row is
independently replayed from the pinned binary.

[COMPASS](https://github.com/rohitsinghlab/compass) separates shared and
intervention-specific response components. Its documented estimators require
measured perturbation responses in the target cell line to fit the shared axis.
Our competition supplies target controls only, so we cannot directly use that
API or silently substitute hidden target responses. A control-only extension
would need its own validation.

## Fixed evidence already obtained

|Mean error ratio to V24, lower better|V38 hurdle25|V39 KOLF hurdle25|V40 rank25|
|---|---:|---:|---:|
|LFC MSE|0.967775|0.967862|0.992838|
|LFC absolute error|0.984076|0.984085|0.998862|
|CPM squared error|0.918643|0.918223|0.980231|

V38's zero-head control has MSE ratio0.966839 and CPM-error ratio0.917319,
slightly better than its learned head.

|Exact development count metric|V24|V38|V39|
|---|---:|---:|---:|
|PDS, higher better|0.791667|0.708333|0.708333|
|Fidelity, higher better|0.581703|0.580967|0.580751|
|Reach, higher better|0.194389|0.169159|0.188606|
|Jaccard, higher better|0.070027|0.069337|0.069581|
|NMAE, lower better|0.912601|0.914119|0.911955|
|Expression error ratio, lower better|37.471936|35.838033|35.829206|

**V38 and V39 fail their frozen count guards and are not released.** HCT116's
expression denominator is very small and its raw error ratio very large; the
table's expression values are not host anchor-scaled scores. HCT116 NMAE is
undefined and is never filled with zero; Orion remains represented by HEK293T.

V40's44-case ideal pooled-profile PDS screen is0.733856 versus0.731782 for
V24 and0.728806 for the contrastive-label-shuffled model. Integer counts and
six-metric confirmation remain separate checks; these ideal margins are not
a private-score estimate.

## KOLF relevance audit

Average per-target pairwise-finite Spearman agreement, no pseudocount:
H1=0.07737, CD4=0.00761, HCT116=0.00828, HEK293T=0.00736,
K562=0.00578. H1 active-response sign agreement is0.6322 versus0.5041 under
target shuffle. Agreement is weaker in the other lineages. Thus broad target
coverage does not imply broad context transfer. V39 tests the source as
auxiliary evidence, not as invented fold changes.

## Validation and operational history

All source/model/input checksums are retained. Tests cover excluded-study
mutation invariance, unknown values, own-target masks, raw integer libraries,
constant-column estimator serialization, every KOLF extracted target,
positive bounded gain, masked-label loss/gradient invariance, and correct
versus wrong identity contrasts.

V38 initially hit a sklearn/NumPy all-missing-column binning error before any
fit. The documented adapter drops columns constant in **training** values
only. Its first untrained launch was stopped to use an importable serialized
class. The original plan and model recipe were preserved.

Transport was withdrawn after severe held-treatment fidelity losses, before
opening additional24 count cases. Original code/plans/results are preserved;
the amended runner changes only the arms evaluated, not any surviving
prediction, seed, split, coefficient or numerical generator. Exact zero-effect
and V24 comparator populations are reused only when the entire count function
and inputs are identical. No model-dependent treatment populations are reused.

The additional24 count cases have previously spent mean/source information.
They are count-level algorithm confirmation, not blind biological validation.
The separate D3 reserve stays unopened.

## Reproduce

Use Python3.12 and the pinned scientific dependencies from
requirements-v35-science.txt. The official competition environment is separate
and is used for packaging.

```cmd
cd /d L:\promec\Animesh\Download\Virtual
set "PY=C:\Users\animeshs\AppData\Local\Temp\vcc_v20_stack\env\Scripts\python.exe"
"%PY%" test_v38_response_hurdle.py
"%PY%" -c "import v38_binning_adapter as a; a.main()"
"%PY%" run_v38_response_hurdle.py prepare
"%PY%" run_v38_response_hurdle.py means
"%PY%" run_v38_response_hurdle.py null
"%PY%" run_v38_standard_counts.py development --version 38
"%PY%" run_v38_response_hurdle.py select
"%PY%" collect_v39_kolf_features.py
"%PY%" test_v39_kolf_features.py
"%PY%" analyze_v39_kolf_support.py
"%PY%" run_v39_kolf_hurdle.py train
"%PY%" run_v39_kolf_hurdle.py prepare
"%PY%" run_v39_kolf_hurdle.py means
"%PY%" reuse_v39_comparators.py
"%PY%" run_v38_standard_counts.py development --version 39
"%PY%" run_v39_kolf_hurdle.py select
"%PY%" test_v40_identity_calibration.py
"%PY%" v40_identity_calibration.py train
"%PY%" v40_identity_calibration.py train --shuffle
"%PY%" run_v40_identity_calibration.py prepare
"%PY%" run_v40_identity_calibration.py means
"%PY%" run_v40_identity_calibration.py null
"%PY%" run_v40_identity_calibration.py development
"%PY%" run_v40_identity_calibration.py select
```

Freeze is a one-time step in each fresh namespace; frozen existing models and
plans refuse overwrite. Set a separate replay output directory for reruns.
Scripts, plans, receipts and small outcomes are versioned; large raw files,
models, count matrices and upload packages stay outside Git. No script uploads
a submission or pushes Git.

