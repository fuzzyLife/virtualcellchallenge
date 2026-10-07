"""Write a reproducible pilot report and scientific comparison figure."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import build_submission as build
from train_embedding_response import ROOT, OUT


def table(frame):
    rows=['| '+' | '.join(frame.columns)+' |','| '+' | '.join(['---']*len(frame.columns))+' |']
    for row in frame.itertuples(index=False,name=None):
        rows.append('| '+' | '.join(f'{v:.6f}' if isinstance(v,(float,np.floating)) else str(v) for v in row)+' |')
    return '\n'.join(rows)


def main():
    decision=json.loads((OUT/'adoption_decision.json').read_text())
    validation=json.loads((OUT/'validation.json').read_text());assert validation['passed'] and validation['count_evaluation_complete']
    contrast=pd.read_csv(OUT/'contrast_summary.csv');cells=pd.read_csv(OUT/'generated_cell_summary.csv')
    source=pd.read_csv(OUT/'source_rows.csv')
    variants=['identity_rank','pca_rank','tahoe_rank','tahoe_pool']
    fig,axes=plt.subplots(1,2,figsize=(12,4.3))
    colors=['#1f6f8b','#d08720','#33804b','#8464a4']
    for j,variant in enumerate(variants):
        group=contrast[(contrast.regime=='lineage_holdout')&(contrast.variant==variant)].set_index('family').reindex(['K562','CD4','H1'])
        axes[0].bar(np.arange(3)+(j-1.5)*.18,100*(group.mse_ratio_to_zero-1),width=.18,label=variant,color=colors[j])
    axes[0].axhline(0,color='black',linestyle='--',linewidth=1)
    axes[0].set_xticks(range(3),['K562','CD4','H1']);axes[0].set_ylim(-2,1)
    axes[0].set_ylabel('Contrast MSE change vs zero-response (%)')
    axes[0].set_title('Entire source lineage held out');axes[0].legend(fontsize=7)
    order=['v4_transfer']+variants;indexed=cells.set_index('candidate').loc[order]
    axes[1].bar(np.arange(5)-.16,indexed.fidelity,width=.32,label='Fidelity',color='#1f6f8b')
    axes[1].bar(np.arange(5)+.16,indexed.jaccard,width=.32,label='DE overlap',color='#d08720')
    axes[1].set_xticks(range(5),['V4 transfer','Rank','PCA','Tahoe/rank','Tahoe/pool'],rotation=25,ha='right')
    axes[1].set_ylabel('Raw source metric (higher is better)');axes[1].legend(fontsize=8)
    axes[1].set_title('H1 generated cells; 25% response blends')
    fig.suptitle('Interventional embedding pilot: public-source diagnostics, not leaderboard scores',fontsize=11)
    fig.tight_layout();fig.savefig(OUT/'pilot_comparison.png',dpi=170);plt.close(fig)
    base=cells[cells.candidate=='v4_transfer'].iloc[0]
    rows=[]
    for row in cells.itertuples(index=False):
        rows.append({'Candidate':row.candidate,'Fidelity':row.fidelity,'DE overlap':row.jaccard,
            'Expression ratio':row.expression_ratio,'PDS':row.pds,
            'Profile MSE':row.uncorrected_profile_mse,'Variance / template':row.variance_ratio_to_template})
    text='''# Interventional embedding response pilot

This pilot implemented the idea recommended in `CROSS_MODAL_MODEL_REVIEW.md`: use frozen biological gene priors plus control-derived context features to predict measured CRISPRi responses. It trained small masked kernel-ridge heads, evaluated whole-lineage and whole-target holdouts, and tested the resulting priors through V4's count generator.

**Decision: retain V4.** None of the tested priors met the adoption criteria. No new upload package was generated. This conclusion concerns these small response heads and cached RNA-only embeddings; it does not rule out a properly trained cross-modal model with more interventional data.

## Inputs and source preparation

The training set contains {nrows} target/context contrasts spanning {ntargets} unique panel targets:

{source_table}

- K562: the existing Replogle raw pseudobulk source; aggregate effective constructs (at least 30 cells, source knockdown at least 50%), sum duplicate symbols, estimate clipped log2 response and Poisson noise. Regression labels subtract 1.96 noise standard errors. Normal-approximation K562 calls are diagnostics, not observed cell-level DE calls.
- CD4: quality-filtered merged donor response summaries, retaining Rest/Stim8hr/Stim48hr. Use source log2 fold changes, standard errors and adjusted p-values. Fixed regression-label shrinkage is `0.04 / (0.04 + SE^2)` as in the previous model. A CD4 holdout excludes all three conditions together.
- H1: the existing 25-target raw-cell contrast cache, using full public response summaries only when H1 is a training lineage. Labels apply the previous `0.7 * sign(LFC) * max(abs(LFC)-SE, 0)` rule. All H1 response rows are removed when testing H1 transfer.

Source hashes, masks, row identities and selected group sizes are recorded. This is a selected public development panel, not an unbiased sample of the challenge distribution. The public data were previously inspected. Some merged source summaries share targets/donors; lineage-level exclusions prevent them crossing the lineage test boundary.

## Models

All four models use the same output masks, family-balanced regression loss and fixed ridge penalty **10**. There was no outcome-based penalty or blend search. Training response means are clipped to +/-3 log2 units. A separate solve for each observation-mask pattern excludes unavailable labels; a missing response is never treated as measured zero. Unavailable response/error values are masked per row. H1 output coverage is selected using non-targeting control CPM >5, without consulting test-target responses.

The kernel combines a same-intervention indicator, a weak global intercept, an optional gene-vector similarity and a control-context similarity. Gene-vector centering/scaling uses distinct training targets only. A radial-basis context kernel modulates intervention similarity; a shared context component permits conservative transfer.

| Variant | Intervention prior | Context features |
|---|---|---|
| `identity_rank` | Known target identity; no pretrained vectors | Pairwise-valid Spearman similarities to training controls |
| `pca_rank` | Cached 15-dimensional DepMap PCA gene vectors | Same rank features |
| `tahoe_rank` | Cached 512-dimensional RNA-only Tahoe gene vectors | Same rank features |
| `tahoe_pool` | Same Tahoe gene vectors | Engineered rank-weighted pooling of gene vectors from controls |

Context rank correlations remove nonfinite pairs before ranking, retain observed RNA zeros and use average ranks for ties. No pseudocount is added to the rank features. Missing gene embeddings explicitly remove that gene's embedding contribution, preserving its identity/global fallback. Engineered pooling is an exploratory feature; it is not the unavailable Trove cell encoder.

All 300 panel genes are excluded as fitted output columns. They retain the baseline response when count blends are evaluated; target identity remains an input feature. This avoids own-target leakage and focuses the pilot on downstream genes. Counts and metrics still use the complete H1 gene axis, with the own target excluded from direction/expression metrics and all panel targets excluded from discrimination.

## Validation design

Six outer splits were evaluated: three entire-lineage holdouts (K562, CD4, H1) and three target-group folds. Every target response is absent across all contexts in a target holdout. Every response from the held-out lineage is absent in a lineage holdout. Unperturbed test-context controls are allowed as query features. No query uses a measured treated expression profile.

The summary contrast evaluation uses identical gene universes for all representations in a split, source-provided standard errors and a normal/BH call approximation. It reports proxy fidelity/overlap, not official scoring. Empty-DE rows can receive overlap 1, so average proxy overlap must be read alongside the call counts. K562, CD4 and H1 have different uncertainty/calibration conventions; their diagnostics are reported separately.

The independent validator checks all 7,112 contrast records, exact hashes, source masks and fold boundaries. During verification, condition-wide missingness selection was replaced with controls-only coverage and per-row masks to protect target folds. The H1-held-lineage prediction cache remained byte-for-byte identical; `mask_boundary_verification.json` records that equality, so its completed cell evaluation was reused. Four focused tests check missing-label exclusion, observed-zero rank semantics, train/query separation and unsupported output handling.

## Whole-lineage contrast results

{contrast_table}

The heads largely predict weak responses. On H1, Tahoe/rank has slightly lower MSE than the rank-only head, but worse sign agreement on the measured DE subset. Neither Tahoe representation improves error and DE diagnostics consistently across all three held-out lineages. Small error changes around the zero-response baseline therefore do not justify adoption.

## Generated-cell evaluation

Evaluated {primary_count} H1 target groups with at least 200 actual cells in both batch partitions. Each comparison uses up to 400 distinct target/template cells, two seeds and two partitions. Non-targeting template and reference controls are disjoint by source cell ID. All five candidates use exactly the same actual cells, template rows, random seed, 90% direct knockdown and V4 zero-aware binomial/Poisson generator.

The baseline is the frozen V4 **cross-lineage transfer** formula `{cd4_gain} * CD4_specific_C + {k562_gain} * K562`. ALL H1 direct and shared-prior components are removed. The CD4 condition mixture is V4's C mixture; it is a fixed development comparator, not a claim that H1 is context C. These transfer coefficients were previously calibrated with H1 development data. The head holdout is strict, but the experiment as a whole is not a fresh blind benchmark. It does not reproduce or estimate the exact published A/B/C V4 artifact.

Each head supplies a fixed 25% replacement of the baseline response on supported non-panel coordinates; unsupported coordinates and panel outputs retain the baseline. Query-control CPM expression gates are applied consistently after blending. Response signs may change and previously zero coordinates may become nonzero: these are tested rather than automatically protected.

Actual generated cells are evaluated with Wilcoxon/BH rank-DE, direction fidelity, reach, DE overlap, log-FC error, panel discrimination, expression correction, detection and dispersion. These are raw, unanchored public-source metrics, not leaderboard-normalized values. The expression ratio is a source-panel capped-jackknife ratio; it has no competition anchor set.

{count_table}

Paired bootstrap intervals resample the {primary_count} target groups (2,000 draws), averaging seeds/partitions within a target. This captures target-to-target variation within this selected panel, not uncertainty across unseen contexts. The adoption gate requires consistent source-summary behavior plus supported fidelity, overlap and expression improvements, preserved discrimination and at least 90% of baseline dispersion. Detailed intervals and every failed guard are in `adoption_decision.json`.

The blends reduce uncorrected profile error but worsen the corrected expression ratio, DE fidelity, overlap and discrimination. For Tahoe/rank, fidelity falls from 0.145431 to 0.057594 and overlap from 0.014226 to 0.008027; their paired bootstrap difference intervals are entirely negative. The expression ratio rises from 1.276690 to 1.283540. Cell variance remains close to the baseline, but the expression correction factor falls from about 0.387 to 0.271 as between-target profile spread shrinks. Less noise correction offsets the small uncorrected-MSE improvement. This repeats the earlier lesson: smoothing an intervention toward a generic weak response can improve average expression error while removing correctly directed significant changes. Retaining cell-template variance alone is insufficient; response strength and sign agreement also matter.

![Pilot comparisons](analysis/embedding_response_pilot/pilot_comparison.png)

## Reproduce

From Windows Command Prompt in the root folder:

```bat
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m unittest test_embedding_response -v
.\\.competition_env\\Scripts\\python.exe train_embedding_response.py
.\\.competition_env\\Scripts\\python.exe evaluate_embedding_response.py
.\\.competition_env\\Scripts\\python.exe validate_embedding_response.py
.\\.competition_env\\Scripts\\python.exe document_embedding_response.py
```

`train_embedding_response.py --reuse-cache` reruns fitting from the verified local cache. All code is in the root. Small CSV/JSON evaluation/provenance files and the figure are tracked; numerical `.npz` caches are ignored and regenerable. No protein checkpoint or training dataset was downloaded, and no competition API submission was made.

The next useful change should address response transfer rather than only reconstruction error: additional context-matched measured perturbations, a response model with stronger interventional conditioning, or a released cross-modal encoder tested with the same controlled generator. This fixed-ridge pilot is a limited model family; lower penalties or other heads remain possible research experiments, not established improvements.
'''
    source_table=source.groupby('context').agg(contrasts=('target_gene','size'),unique_targets=('target_gene','nunique')).reset_index()
    selected=contrast[contrast.regime=='lineage_holdout'][['family','variant','mse_ratio_to_zero','sign_accuracy_reference_DE','fidelity_proxy','jaccard_proxy']].copy()
    text=text.format(nrows=len(source),ntargets=source.target_gene.nunique(),source_table=table(source_table),
        contrast_table=table(selected),primary_count=len(decision['primary_targets']),cd4_gain=.5,k562_gain=1.,count_table=table(pd.DataFrame(rows)))
    (ROOT/'EMBEDDING_RESPONSE_PILOT.md').write_text(text,encoding='utf-8')
    manifest={'report_sha256':build.sha256(ROOT/'EMBEDDING_RESPONSE_PILOT.md'),
        'plot_sha256':build.sha256(OUT/'pilot_comparison.png'),'decision_sha256':build.sha256(OUT/'adoption_decision.json'),
        'documentation_script_sha256':build.sha256(Path(__file__))}
    build.write_json(OUT/'documentation_manifest.json',manifest)
    print('Written EMBEDDING_RESPONSE_PILOT.md; adoption:',decision['selected'])


if __name__=='__main__':main()
