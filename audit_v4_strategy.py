"""Plot V4 source tradeoffs and audit selected generated downstream responses."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from audit_prediction import ControlRankTest,normalized_counts
import build_submission as build

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'experiments/v4_lineage_response';A=OUT/'analysis'

def main():
    model=json.loads((A/'reference_model.json').read_text());manifest,genes,targets=build.read_spec(ROOT)
    with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as z:effects=z['effects'].copy()
    with np.load(ROOT/'experiments/v3_context_ensemble/analysis/reference_model_arrays.npz',allow_pickle=False) as z:old=z['effects'].copy()
    rows=[]
    for ci,context in enumerate(['A','B','C']):
        for ti,target in enumerate(targets):
            rows.append({'context':context,'target_gene':target,'v3_effects_above_0_1':int((np.abs(old[ci,ti])>.1).sum()),
                'v4_effects_above_0_1':int((np.abs(effects[ci,ti])>.1).sum()),'v3_effect_rms':float(np.sqrt((old[ci,ti]**2).mean())),
                'v4_effect_rms':float(np.sqrt((effects[ci,ti]**2).mean()))})
    comparison=pd.DataFrame(rows);comparison.to_csv(A/'v3_v4_effect_comparison.csv',index=False)
    generated=[]
    names=['AKT2','SMARCA5','VIM','S100A11','FXYD5','KIF21B']
    previous=ad.read_h5ad(ROOT/'experiments/v3_context_ensemble/prediction.h5ad',backed='r')
    prediction=ad.read_h5ad(OUT/'prediction.h5ad',backed='r')
    for ci,context in enumerate(['A','B','C']):
        original=ad.read_h5ad(ROOT/f'context_{context}.h5ad');controls=normalized_counts(original.X)
        cmean=np.asarray(controls.mean(axis=0)).ravel();keep=cmean>5;columns=np.flatnonzero(keep)
        rank_test=ControlRankTest(controls[:,columns]);labels=genes[columns]
        for target in names:
            ti=list(targets).index(target);start=(ci*300+ti)*400
            for version,obj in [('v3',previous),('v4',prediction)]:
                block=obj.X[start:start+400].tocsr();normalized=normalized_counts(block)[:,columns]
                p,q,lfc=rank_test.test(normalized);downstream=labels!=target;sig=(q<.05)&downstream
                generated.append({'context':context,'target_gene':target,'version':version,'downstream_DE_calls':int(sig.sum()),
                    'up_calls':int((sig&(lfc>0)).sum()),'down_calls':int((sig&(lfc<0)).sum()),
                    'mean_library':float(block.sum()/400),'cell_count':400})
        print(context,'selected generated response comparisons audited',flush=True)
        original.file.close();del original,controls,rank_test
    previous.file.close();prediction.file.close();actual=pd.DataFrame(generated);actual.to_csv(A/'v3_v4_selected_response_audit.csv',index=False)
    build.write_json(A/'response_audit_scope.json',{'targets':names,'groups_each_version':18,'method':'Original controls; CPM >5 filter; cached exact two-sided Wilcoxon with ties, continuity and within-target BH; own target excluded',
        'limitation':'Counts of predicted response genes measure model activity, not correct directions, recall, overlap or leaderboard scores; provided controls differ from hidden scoring controls',
        'script_sha256':build.sha256(Path(__file__))})
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    cd4=model['cd4_selection'];h1=model['h1_selection'];cross=model['cross_lineage_selection']
    x=np.arange(3);labels=['CD4 donors','H1 direct','Cross-lineage H1']
    bases=[cd4['v3_fixed_baseline'],h1['v3_fixed_baseline'],cross['zero_baseline']]
    vals=[cd4['nested_procedure'],h1['nested_procedure'],cross['nested_procedure']]
    axes[0,0].bar(x-.18,[v['direction_fidelity_proxy'] for v in bases],width=.36,label='Comparison baseline')
    axes[0,0].bar(x+.18,[v['direction_fidelity_proxy'] for v in vals],width=.36,label='Nested selection procedure')
    axes[0,0].set_xticks(x,labels);axes[0,0].set_ylabel('Direction-yield proxy');axes[0,0].set_title('Independent source contrasts; 400-cell normal call proxy');axes[0,0].legend(fontsize=8)
    axes[0,1].bar(x,[v['relative_zero_loss'] for v in vals],color=['#4477aa','#228833','#ee7733'])
    axes[0,1].axhline(1,c='black',ls='--',label='Zero downstream effect');axes[0,1].set_ylim(.75,1.08)
    axes[0,1].set_xticks(x,labels);axes[0,1].set_ylabel('Log2FC squared error / zero-effect error');axes[0,1].set_title('Coverage comes with an explicit error tradeoff');axes[0,1].legend(fontsize=8)
    sums=comparison.groupby('context')[['v3_effects_above_0_1','v4_effects_above_0_1']].sum().reindex(['A','B','C'])
    axes[1,0].bar(x-.18,sums.iloc[:,0],width=.36,label='V3');axes[1,0].bar(x+.18,sums.iloc[:,1],width=.36,label='V4')
    axes[1,0].set_xticks(x,['A','B','C']);axes[1,0].set_ylabel('Modeled gene-target effects');axes[1,0].set_title('Effects exceeding |log2FC| 0.1; not measured DE');axes[1,0].legend()
    parts=actual.groupby(['context','version']).downstream_DE_calls.mean().unstack().reindex(['A','B','C'])
    axes[1,1].bar(x-.18,parts.v3,width=.36,label='V3');axes[1,1].bar(x+.18,parts.v4,width=.36,label='V4')
    axes[1,1].set_xticks(x,['A','B','C']);axes[1,1].set_ylabel('Mean downstream significant calls');axes[1,1].set_title('Six selected targets against provided controls; no accuracy claim');axes[1,1].legend()
    fig.suptitle('V4: lineage-guided responses, coverage calibration and independent control draws',fontsize=14)
    fig.savefig(A/'v4_strategy_audit.png',dpi=180);fig.savefig(A/'v4_strategy_audit.pdf');plt.close(fig)
    print(actual.to_string(index=False),flush=True)
    document_results()


def document_results():
    ready=json.loads((A/'submission_ready.json').read_text())
    independent=json.loads((A/'independent_validation.json').read_text())
    model_validation=json.loads((A/'model_validation.json').read_text())
    model=json.loads((A/'reference_model.json').read_text())
    assert ready['official_prep_passed'] and independent['passed'] and model_validation['passed']
    # Readiness preserves generation-time status; publication is recorded separately.
    publication_path=A/'leaderboard_result.json'
    publication=json.loads(publication_path.read_text()) if publication_path.exists() else None
    if publication:
        assert publication['final_status']=='published' and publication['package_sha256']==ready['sha256']
    status_text=(f"Published as `{publication['entry_id']}` with overall score **{publication['overall']:.6f}**. "
        f"Rank was {publication['rank_at_publication']} at publication and {publication['rank_at_api_check']} at the recorded API check; ranks change over time."
        if publication else 'It has not been uploaded or scored.')
    records=pd.read_csv(A/'target_response_sources.csv')
    response=pd.read_csv(A/'v3_v4_selected_response_audit.csv')
    text=('# V4: published lineage-guided response submission\n\n' if publication else '# V4: upload-ready lineage-guided response submission\n\n')
    text+='The new candidate is generated, officially packaged and independently validated. '+status_text+' The package is isolated in `experiments/v4_lineage_response/`; the earlier submission files retain their original locations.\n\n'
    text+='| Check | Completed result |\n|---|---|\n'
    items=[('Package','`experiments/v4_lineage_response/prediction.vcc`'),('Package bytes',f"{ready['bytes']:,} ({ready['bytes']/1024**3:.2f} GiB)"),
        ('H5AD bytes',f"{(OUT/'prediction.h5ad').stat().st_size:,}"),('Shape','360,000 cells x 18,533 genes'),('Groups','900, exactly 400 cells each'),
        ('Stored count entries',f"{ready['prediction_validation']['nnz']:,}"),('Maximum cell total',f"{ready['prediction_validation']['maximum_cell_counts']:,}"),
        ('Official CLI preparation',f"Passed; CLI {ready['package_meta']['cli_version']}"),('Regression tests','23 passed'),
        ('Original source control rows','All 360,000 verified; 400 unique rows per group; 46 balanced guides'),
        ('Independent count expectation checks',f"{independent['independent_gene_count_expectations']:,} across {len(independent['independent_count_expectation_groups'])} groups, with complete deterministic count replay"),
        ('All final model effects','16,679,700 independently assembled'),('Generation and packaging time',f"{ready['elapsed_seconds']/60:.1f} minutes"),('Leaderboard result',f"Published; overall {publication['overall']:.6f}" if publication else 'Not yet available')]
    text+=''.join(f'| {key} | {value} |\n' for key,value in items)
    text+='\nPackage SHA-256:\n\n```text\n'+ready['sha256']+'\n```\n\n'
    published_section=''
    if publication:
        comparison=pd.read_csv(A/'leaderboard_comparison.csv')
        published_section='## Published validation result\n\n'
        published_section+='The first, V2 and V4 scores use the same validation panel and anchor set. No V3 submission score is recorded, so there is no V3 leaderboard comparison.\n\n'
        published_section+='| Scaled component | First | V2 | V4 | Change from V2 |\n|---|---:|---:|---:|---:|\n'
        for row in comparison.itertuples():
            published_section+=f'| {row.metric} | {row.first:.6f} | {row.v2:.6f} | {row.v4:.6f} | {row.change_from_v2:+.6f} |\n'
        raw=publication['raw_metrics'];previous=json.loads((ROOT/'experiments/v2_cd4_hybrid/analysis/published_api_status.json').read_text())
        published_section+=f"\nOverall improved by **{publication['overall_change_from_v2']:+.6f}** versus V2 and **{publication['overall_change_from_first']:+.6f}** versus the first entry. All five non-expression scaled components improve. Direction fidelity accounts arithmetically for {publication['fidelity_fraction_of_overall_gain_from_v2']:.1%} of the overall gain; this is a score contribution, not causal attribution to a model change.\n\n"
        published_section+=f"Raw expression error worsens from {previous['expr_mse_unbiased_capped_norm']:.6f} to {raw['expr_mse_unbiased_capped_norm']:.6f} ({publication['raw_expression_error_relative_change_from_v2']:+.1%}), although the scaled expression score stays at zero. Zero is a clamped leaderboard score, not zero prediction error. Raw direction-fidelity yield rises from {previous['de_wilcoxon_direction_fidelity_yield_raw']:.6f} to {raw['de_wilcoxon_direction_fidelity_yield_raw']:.6f}; it is response-budget-normalized yield rather than sign precision. Raw significant-set Jaccard rises from {previous['de_wilcoxon_sig_jaccard']:.6f} to {raw['de_wilcoxon_sig_jaccard']:.6f}, leaving substantial response-set mismatch.\n\n"
        published_section+='V4 is the best of the three published project submissions on this validation panel. Overall remains below the zero-score anchor, with direction fidelity the largest negative component. The package changes lineage policy, effect strength, transfer priors and control sampling together; aggregate scores cannot isolate which change helped, which context remains weak, or whether residual fidelity loss comes from missed responses or wrong directions. The next experiments should balance response coverage with expression error using source validation and controlled ablations.\n\n'
        published_section+='Records: `leaderboard_result.json`, `published_api_status.json`, `leaderboard_comparison.csv` and `raw_metric_comparison.csv` in the candidate analysis folder. `submission_ready.json` preserves the generation-time snapshot; the separate publication record records completed upload and scoring. Upload took 13 minutes 20 seconds. [Official scoring definitions](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics.md).\n\n'
        published_section+='Record or refresh the read-only API confirmation and regenerate these tables with the root script `record_v4_result.py`:\n\n```cmd\n.\\.competition_env\\Scripts\\python.exe record_v4_result.py\n```\n\nThe supplied publication rank is preserved separately from the current API rank. This command does not upload, regenerate predictions, or alter the package.\n\n'
        text+=published_section
    text+='## Strategy\n\nA uses supported CD4 responses with a donor-calibrated target-specific gain of 1 instead of V3\'s 0.5. K562 fills unsupported coordinates. B/C use measured H1 responses for 13 targets; other targets use 0.5 CD4 target-specific effect + the existing 0.3-shrunk K562 effect + 0.5 shared H1 response. These are response-estimator coefficients, not cell-type fractions. Gene/target expression compatibility moderates weakly expressed coordinates. Every group starts from a separately seeded, guide/depth-stratified draw of 400 original controls.\n\n'
    text+='Nested CD4 selection increases the direction-yield proxy from 0.02324 to 0.17701, with predicted response calls increasing from 2.81 to 31.86 per source contrast. Relative squared error changes from 0.98425 to 1.00436. Public H1 direct validation retains the V3 estimator (relative squared error 0.80299). The cross-lineage selection procedure has direction-yield proxy 0.40430 and relative squared error 1.02142; the final fixed candidate has 1.01191. These are source-domain tradeoffs using a stated 400-cell normal call proxy, not competition scores or proof of improved B/C predictions.\n\n'
    text+='All 284 previously covered targets retain modeled downstream effects in A. B/C have response estimates for all 300 through the shared prior, including targets with no direct source measurements. The H1 prior is based on only 13 selected interventions and may not represent the hidden panel. More predicted response genes do not establish correct signs or a better leaderboard result.\n\n'
    text+='Complete derivation, objective, fold separation, source limitations and reproduction: [V4_MODEL_STRATEGY.md](V4_MODEL_STRATEGY.md).\n\n'
    text+=('## Completed upload command (Windows Command Prompt)' if publication else '## Submit from Windows Command Prompt')+'\n\n'+('The command below already succeeded; rerunning it creates another submission.' if publication else 'Reuse the existing VCC login:')+'\n\n```cmd\ncd /d L:\\promec\\Animesh\\Download\\Virtual\n.\\.competition_env\\Scripts\\python.exe -m vcc submit .\\experiments\\v4_lineage_response\\prediction.vcc -m "Lineage-guided response coverage v4" --wait\n```\n\n'
    text+='## Actual generated-response diagnostic\n\nThe table compares six deliberately selected targets in each context against the provided controls using the control-mean >5 CPM filter, two-sided Wilcoxon and within-target BH, excluding the own target. It measures predicted response activity rather than hidden response accuracy.\n\n'
    from analyze_target_rank_correlations import table
    text+=table(response)+'\n\n![V4 strategy and diagnostics](experiments/v4_lineage_response/analysis/v4_strategy_audit.png)\n\n'
    text+='## Code and records\n\nRoot scripts are `fit_lineage_submission.py`, `test_lineage_submission.py`, `validate_v4.py`, `audit_v4_strategy.py`, with generation integration in `build_submission.py`. The audit script also regenerates this results document and the README section after package validation.\n\n'
    text+='Tracked records in the candidate\'s `analysis/` folder include nested source-selection JSON/CSV tables, source hashes, response coverage, control/gene summaries, all 360,000 source-row records, model and package validation, exact count-replay evidence, official preparation output, run settings, V3/V4 response diagnostics and PNG/PDF figures. Large source-component caches, effect caches, H5AD and VCC files are excluded from Git and rebuilt locally. No new public source download or leaderboard submission occurs in these scripts.\n'
    (ROOT/'V4_SUBMISSION.md').write_text(text,encoding='utf-8')
    readme=ROOT/'README.md';text=readme.read_text(encoding='utf-8');heading='\n## V4 lineage-guided response coverage candidate\n'
    if heading in text:text=text[:text.index(heading)]
    text+=heading+'\nThe completed candidate is ready at `experiments/v4_lineage_response/prediction.vcc` ('+f"{ready['bytes']:,} bytes; {ready['bytes']/1024**3:.2f} GiB"+'). Official preparation and independent package verification passed for all 900 groups, 360,000 cells and the exact 18,533-gene axis. '+status_text+' [V4_SUBMISSION.md](V4_SUBMISSION.md) records the package hash, completed checks, actual generated-response audit and Command Prompt upload command.\n\n'
    text+='The model uses the marker-supported T-cell-like program to prefer CD4 responses for A, strengthens donor-supported target effects to address response undercoverage, and uses public H1 direct responses plus held-target-calibrated transfer/shared responses for B/C. It applies explicit basal gene/target expression gates and independently draws 400 real controls per target with guide/depth balance. Source selection improves direction-coverage proxies with a documented squared-error tradeoff; it does not establish a better hidden leaderboard score. [V4_MODEL_STRATEGY.md](V4_MODEL_STRATEGY.md) documents every step and exact reproduction.\n\n'
    text+='New root code is `fit_lineage_submission.py`, `test_lineage_submission.py`, `validate_v4.py` and `audit_v4_strategy.py`; `build_submission.py` adds the isolated `lineage-response` mode and complete original-row provenance. All 23 tests passed. Independent validation recalculates all public-source response components, reassembles all 16,679,700 model effects, verifies all 360,000 original control-row/guide/depth records and checks '+f"{independent['independent_gene_count_expectations']:,}"+' gene-count expectations with exact selected-cell generation replay. Model/package provenance and evidence are retained in the candidate analysis folder.\n'
    if publication:
        text+='\n'+published_section
        old='Two submissions have been uploaded, scored, and published, both with weak predictive performance. A third candidate now uses a context-specific source ensemble with independent source validation; see [V3_SUBMISSION.md](V3_SUBMISSION.md).'
        new=f"Three submissions (first, V2 and V4) have been published. V4 is the best so far on the validation panel, with overall score {publication['overall']:.6f}; direction fidelity and expression error remain weaknesses. No V3 submission score is recorded. See [V4_SUBMISSION.md](V4_SUBMISSION.md) for the published comparison and reproducibility."
        text=text.replace(old,new)
    readme.write_text(text.rstrip()+'\n',encoding='utf-8')

if __name__=='__main__':main()
