"""Write V7 reproducibility/report artifacts only after independent package verification."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import build_submission as build
ROOT=build.ROOT;A=ROOT/'experiments/v7_multisource_adaptive/analysis'

def table(frame):
    def show(v):
        if isinstance(v,(float,np.floating)):return f'{v:.6f}' if np.isfinite(v) else 'n/a'
        return str(v).replace('|','/')
    return '\n'.join(['| '+' | '.join(map(str,frame.columns))+' |','| '+' | '.join(['---']*len(frame.columns))+' |',*['| '+' | '.join(show(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None)]])

def main():
    ready=json.loads((A/'submission_ready.json').read_text());verified=json.loads((A/'independent_validation.json').read_text());assert verified['passed'] and ready['sha256']==verified['package_sha256']
    model=json.loads((A/'reference_model.json').read_text());cell=json.loads((A/'cell_selection.json').read_text());source=json.loads((A/'source_provenance.json').read_text())
    frame=pd.read_csv(A/'cell_summary.csv');confirmation=frame[frame.split=='confirmation'];means=confirmation.groupby(['family','candidate']).mean(numeric_only=True)
    balanced=means.groupby('candidate').mean().reset_index()
    columns=['candidate','pds','fidelity','reach','jaccard','nmae','expression_ratio','template_variance_ratio']
    rows=[];audit={}
    for line in ['HCT116','HEK293T']:
        directory=ROOT/f'reference/v7_sources/{line}';record=json.loads((directory/'extraction.json').read_text());gap=json.loads((directory/'expression_fragment_audit.json').read_text())
        rows.append({'source':line,'selected_cells':record['cells'],'targets':record['targets'],'eligible_targets':int(source['row_counts'][line]),'batches':record['batches'],'measured_genes':record['measured_genes'],'excluded_missing_cells':record['selection_record']['missing_selected_cells']})
        audit[line]={'extraction':record,'missing_cell_ranges':gap['missing_cell_ranges'],'missing_expression_records':gap['missing_expression_records'],'audit_sha256':build.sha256(directory/'expression_fragment_audit.json')}
    pd.DataFrame(rows).to_csv(A/'source_coverage.csv',index=False)
    build.write_json(ROOT/'analysis/v7_source_inventory/completed_source_audit.json',audit)
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    metrics=[('fidelity','Fidelity (higher is better)'),('jaccard','DE overlap (higher is better)'),('pds','Discrimination (higher is better)'),('expression_ratio','Corrected expression ratio (lower is better)')]
    names=balanced.candidate.tolist();colors=['#24749c' if n==model['candidate'] else '#b75a45' if n=='v4_transfer' else '#969da5' for n in names]
    for ax,(metric,label) in zip(axes.ravel(),metrics):
        ax.bar(np.arange(len(names)),balanced[metric],color=colors);ax.set_title(label);ax.set_xticks(np.arange(len(names)),[n.replace('_','\n') for n in names],fontsize=8);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
    fig.suptitle('V7 public-source confirmation comparison; equal weight to three lineages',fontsize=13)
    fig.savefig(A/'source_confirmation_comparison.png',dpi=160);plt.close(fig)
    details=[]
    for line in ['H1','HCT116','HEK293T']:
        scores=pd.read_csv(A/f'cell_scores_{line}.csv');unique=scores.drop_duplicates(['target_gene','seed']);details.append({'source':line,'targets':unique.target_gene.nunique(),'calibration_targets':unique[unique.split=='calibration'].target_gene.nunique(),'confirmation_targets':unique[unique.split=='confirmation'].target_gene.nunique(),'observed_cells_min':int(unique.observed_cells.min()),'observed_cells_max':int(unique.observed_cells.max()),'prediction_cells':400})
    summary=table(balanced[columns]);coverage=table(pd.DataFrame(rows));scope=table(pd.DataFrame(details));decision=table(pd.read_csv(A/'cell_selection_scores.csv'))
    text=f'''# V7 submission: expanded measured-source response selection

The upload-ready file is `experiments/v7_multisource_adaptive/prediction.vcc` ({ready['bytes']:,} bytes). It has **360,000 cells x 18,533 genes**, with all 300 targets x 400 cells in A, B and C. Official packaging and independent full-payload validation passed. The package has not been uploaded and has no hidden leaderboard score.

SHA-256: `{ready['sha256']}`.

## Frozen decision

Generated-cell candidate: **{model['candidate']}**. Count generator: **{model['generator']}**. Production pooling configuration: `{model['production_summary_config']['name']}`. Calibration winner: `{cell['calibration_winner']}`; confirmation improvement accepted: `{cell['confirmed_improvement']}`.

{model['A_policy']}

This uses public CRISPRi interventions, conservative measured-label masks and each query context's own controls. Count sampling seed is 20261008. The 90% own-target knockdown and low-expression gates remain assumptions. The model does not infer the challenge's hidden perturbation truth from control correlations.

## Source coverage and integrity

{coverage}

Both converted-source expression fragments were audited for absent cells and shifted offsets. Duplicate-name suffixes were mapped through the original canonical integer gene metadata and summed, rather than losing counts. Every selected raw cell ID and full library was checked before projection. Two independent original-Parquet checks verified all challenge-axis counts for 32 Batch1 cells per source. K562's full raw matrix was inspected and confirmed to share its existing 7,681-gene challenge overlap; downloading it would not broaden the output axis.

## Public generated-cell evidence

{scope}

There are two paired seeds. Each evaluated lineage is excluded from response fitting; tuning is nested in other lineages. Actual treated cells number 200-400, while all candidates generate 400 predictions and use 400 distinct matched-batch reference controls. Templates/reference controls never overlap. Calibration and confirmation targets are distinct. Previously inspected public data remain development data.

The following averages use only confirmation targets and give equal weight to the three source lineages. Values are raw metrics, **not leaderboard-normalized scores**. A local improvement does not guarantee a hidden score gain.

{summary}

![Public confirmation comparison](experiments/v7_multisource_adaptive/analysis/source_confirmation_comparison.png)

Selection and guards:

{decision}

The utility is a declared mean of raw discrimination, fidelity, reach, overlap, clipped 1-NMAE and clipped 1-expression ratio. Fidelity must remain within 0.01 of baseline, expression ratio within 1.05 times baseline and template variance at least 0.9 times baseline. No candidate is accepted by expression MSE alone. Fast summary-stage normal/SE proxies are followed by actual generated-cell Wilcoxon/BH tests.

## Validation and reproducibility

Nine focused scientific tests passed, including measured-zero versus missing-source handling, matched-batch neutrality, expected-count preservation and conservative ridge label masks. Neutral H1 controls produced zero false DE calls for each tested count generator. Full submission auditing checked all 900 groups, all 360,000 cell libraries and every packaged nonzero value/index/pointer; twelve groups were exactly replayed from their recorded controls and seeds. The independent audit is `experiments/v7_multisource_adaptive/analysis/independent_validation.json`.

The complete steps, pinned sources, assumptions and model scope are in [V7_MODEL_STRATEGY.md](V7_MODEL_STRATEGY.md). Root scripts include preparation, source-offset repair, three neural architectures, four ridge penalties, source-held response/count comparisons, conservative coverage filling, production fitting, sparse generation and independent validation. The small neural response heads implement an intervention-conditioned transition idea; this is not a faithful Arc STATE foundation-model reproduction.

## Upload using Windows CMD

```bat
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc whoami
.\\.competition_env\\Scripts\\python.exe -m vcc submit .\\experiments\\v7_multisource_adaptive\\prediction.vcc -m "Expanded measured-source response v7" --wait
```

If waiting is interrupted after upload, resume status polling for the returned entry ID rather than submitting again:

```bat
.\\.competition_env\\Scripts\\python.exe -m vcc status ENTRY_ID --wait
```

The scripts and reproducibility records are committed locally. Large data caches, model arrays and submission packages are ignored by Git and remain in the project folder.
'''
    (ROOT/'V7_SUBMISSION.md').write_text(text,encoding='utf-8')
    readme=ROOT/'README.md';s=readme.read_text(encoding='utf-8-sig');old='V7 is being built from expanded, stratified HCT116/HEK293T intervention coverage, adaptive response pooling and three perturbation-conditioned neural heads. A converted-source gap was found and is being repaired with full cell-ID/library checks. Root scripts, model selection, generated-cell comparisons and packaging are described in [V7_MODEL_STRATEGY.md](V7_MODEL_STRATEGY.md). No V7 package or hidden-score improvement is claimed until its independent audit passes.'
    new=f'V7 is ready for upload at `experiments/v7_multisource_adaptive/prediction.vcc`, using **{model["candidate"]}** with **{model["generator"]}** count emission. Expanded HCT116/HEK293T samples cover all 300 targets; source gaps and duplicate symbols were repaired and checked against original counts. Official packaging and a full independent payload audit passed. V7 has not been uploaded; V4 remains the best published score. See [V7_SUBMISSION.md](V7_SUBMISSION.md) for the public comparison and upload command, and [V7_MODEL_STRATEGY.md](V7_MODEL_STRATEGY.md) for reproduction.'
    assert old in s;s=s.replace(old,new);readme.write_text(s,encoding='utf-8')
    print('V7 reports and README updated from verified artifacts',flush=True)
if __name__=='__main__':main()
