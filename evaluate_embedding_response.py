"""Evaluate response priors using V4 count generation on held-lineage H1 cells.

Fixed 25% blends; identical template/control cells and two fixed seeds. The V4
comparator uses its CD4/K562 cross-lineage transfer policy with ALL H1 response
components removed. It is a source-transfer comparator, not the published A/B/C
V4 score or the earlier opposite-batch oracle estimator.
"""
from pathlib import Path
import json
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as build
from embedding_response_methods import VARIANTS
from train_embedding_response import OUT, ROOT
from v6_methods import moments, de, direction_metrics, profile_jackknife, summarize_profiles

SEED = 20261007
CANDIDATES = ('v4_transfer', 'identity_rank', 'pca_rank', 'tahoe_rank', 'tahoe_pool')


def paired_bootstrap(frame, candidate, key, iterations=2000):
    average=frame.groupby(['target_gene','candidate'])[key].mean().unstack()
    delta=(average[candidate]-average['v4_transfer']).dropna().to_numpy()
    rng=np.random.default_rng(SEED)
    draws=delta[rng.integers(len(delta),size=(iterations,len(delta)))].mean(axis=1)
    return {'difference':float(delta.mean()),'ci_low':float(np.quantile(draws,.025)),
            'ci_high':float(np.quantile(draws,.975)), 'targets':len(delta)}


def main():
    fit=json.loads((OUT/'fit_manifest.json').read_text())
    cache=OUT/'h1_lineage_predictions.npz';assert build.sha256(cache)==fit['prediction_cache_sha256']
    with np.load(cache,allow_pickle=False) as z:pred={k:z[k].copy() for k in z.files}
    source=json.loads((OUT/'source_provenance.json').read_text());assert build.sha256(OUT/'source_cache.npz')==source['cache_sha256']
    with np.load(OUT/'source_cache.npz',allow_pickle=False) as z:
        challenge_genes=z['genes'];challenge_targets=z['targets'];row_targets=z['row_targets'];row_contexts=z['row_contexts']
    directory=ROOT/'reference/v5_h1';extraction=json.loads((directory/'extraction.json').read_text())
    for name,key in [('selected_counts.npz','counts_sha256'),('selected_cells.csv','selected_cells_sha256'),('gene_names.csv','gene_names_sha256')]:
        assert build.sha256(directory/name)==extraction[key]
    x=sparse.load_npz(directory/'selected_counts.npz').astype(np.int32)
    meta=pd.read_csv(directory/'selected_cells.csv');genes=pd.read_csv(directory/'gene_names.csv').gene_name.to_numpy(dtype=str)
    assert meta.source_cell_id.is_unique
    pos={g:i for i,g in enumerate(challenge_genes)};si=np.array([i for i,g in enumerate(genes) if g in pos]);di=np.array([pos[genes[i]] for i in si])
    v4=ROOT/'experiments/v4_lineage_response/analysis'
    record=json.loads((v4/'reference_model.json').read_text());components=v4/'source_components.npz'
    assert build.sha256(components)==record['source_components_sha256']
    with np.load(components,allow_pickle=False) as z:
        np.testing.assert_array_equal(z['genes'],challenge_genes);np.testing.assert_array_equal(z['targets'],challenge_targets)
        transfer=record['cross_lineage_settings']
        baseline=transfer['cd4_coefficient']*z['cd4_specific'][2]+transfer['k562_coefficient']*z['k562']
    effects={name:{} for name in CANDIDATES};embedding_support={}
    for variant in VARIANTS:
        rows=pred[variant+'_rows'];support=pred[variant+'_supported']
        assert (row_contexts[rows]==4).all()
        embedding_support[variant]=int(support.sum())
        for j,row in enumerate(rows):
            ti=int(row_targets[row]);effect=baseline[ti].copy()
            effect[support]=.75*effect[support]+.25*pred[variant][j,support]
            aligned=np.zeros(len(genes),np.float32);aligned[si]=effect[di];effects[variant][str(challenge_targets[ti])]=aligned
    h1_names=sorted(effects['identity_rank'])
    for target in h1_names:
        ti=int(np.flatnonzero(challenge_targets==target)[0]);aligned=np.zeros(len(genes),np.float32);aligned[si]=baseline[ti,di]
        effects['v4_transfer'][target]=aligned
    batches=sorted(meta.batch.unique());partition=meta.batch.map({b:i%2 for i,b in enumerate(batches)}).to_numpy()
    controls=meta.target_gene.to_numpy()=='non-targeting';labels=meta.batch.to_numpy();names=meta.target_gene.to_numpy()
    counts=pd.crosstab(names,partition).reindex(h1_names)
    primary=counts.index[(counts.min(axis=1)>=200)].to_numpy(dtype=str)
    assert len(primary)>=5
    counts.to_csv(OUT/'h1_group_sizes.csv')
    records=[];profiles={};identity=[]
    for part in [0,1]:
        templates=[];references=[]
        for batch in batches:
            rows=np.flatnonzero(controls&(labels==batch)&(partition==part))
            if not len(rows):continue
            rng=np.random.default_rng(SEED+part+sum(str(batch).encode()));rows=rng.permutation(rows)
            templates.extend(rows[:len(rows)//2]);references.extend(rows[len(rows)//2:])
        templates=np.array(templates);references=np.array(references)
        assert not set(templates)&set(references)
        template=x[templates];ref=x[references]
        depth,_,_,tv=moments(template);mean_raw=np.asarray(template.mean(axis=0)).ravel()
        rp,rj=profile_jackknife(ref)
        def evaluate_target(target):
            actual_rows=np.flatnonzero((names==target)&(partition==part));own=int(np.flatnonzero(genes==target)[0])
            ti=int(np.flatnonzero(challenge_targets==target)[0]);entries=[];panels={};ids=[]
            # Gate using only the query context's non-targeting controls.
            control_mean=moments(template)[1]*1e6
            gate=np.clip(control_mean/5,0,1)*np.clip(control_mean[own]/5,0,1)
            for seed_index in [0,1]:
                seed=SEED+part*10000+ti*100+seed_index;rng=np.random.default_rng(seed)
                n=min(400,len(actual_rows),len(templates));ar=rng.choice(actual_rows,n,replace=False);rows=rng.choice(len(templates),n,replace=False)
                actual=x[ar];truth=de(actual,ref,genes);ap,aj=profile_jackknife(actual)
                for candidate in CANDIDATES:
                    effect=effects[candidate][target]*gate
                    factor=np.exp2(np.clip(effect,-3,3));factor[own]=.1
                    generated=build.simulate_counts_zero_aware(template,rows,factor,np.random.default_rng(seed),mean_raw,depth.mean())
                    prediction=de(generated,ref,genes);metrics=direction_metrics(prediction,truth,genes,own)
                    gp,gj=profile_jackknife(generated);_,_,_,gv=moments(generated)
                    valid=np.arange(len(genes))!=own
                    entries.append({'candidate':candidate,'target_gene':target,'partition':part,'seed':seed_index,'cells':n,
                        **metrics,'uncorrected_profile_mse':float(((gp-ap)[valid]**2).mean()),
                        'variance_ratio_to_template':float(gv.sum()/tv.sum()),'mean_detected_genes':float(generated.getnnz(axis=1).mean())})
                    panels[(target,part,seed_index,candidate)]=(gp,ap,rp,gj,aj,rj)
                ids.append({'target_gene':target,'partition':part,'seed':seed_index,'n':n,
                    'actual_source_ids':meta.iloc[ar].source_cell_id.tolist(),
                    'template_source_ids':meta.iloc[templates[rows]].source_cell_id.tolist(),
                    'reference_control_ids':meta.iloc[references].source_cell_id.tolist()})
            return entries,panels,ids
        with ThreadPoolExecutor(max_workers=2) as pool:
            for i,(entries,panels,ids) in enumerate(pool.map(evaluate_target,primary)):
                records.extend(entries);profiles.update(panels);identity.extend(ids)
                pd.DataFrame(records).to_csv(OUT/'generated_cell_scores.csv',index=False)
                print('H1 count evaluation part',part,i+1,'/',len(primary),'targets complete',flush=True)
    frame=pd.DataFrame(records);summary=[];intervals={}
    for candidate in CANDIDATES:
        metric=frame[frame.candidate==candidate].select_dtypes('number').mean().to_dict();panels=[]
        for part in [0,1]:
            for seed in [0,1]:
                arr=[np.stack([profiles[(str(t),part,seed,candidate)][k] for t in primary]) for k in range(6)]
                panels.append(summarize_profiles(*arr,genes,primary,challenge_targets))
        summary.append({'candidate':candidate,**metric,**{k:float(np.mean([p[k] for p in panels])) for k in panels[0]}})
        if candidate!='v4_transfer':
            intervals[candidate]={key:paired_bootstrap(frame,candidate,key) for key in ['fidelity','jaccard','uncorrected_profile_mse','variance_ratio_to_template']}
    summary=pd.DataFrame(summary);summary.to_csv(OUT/'generated_cell_summary.csv',index=False)
    screen=pd.read_csv(OUT/'contrast_summary.csv');base=summary[summary.candidate=='v4_transfer'].iloc[0]
    decisions={}
    for candidate in CANDIDATES[1:]:
        model=summary[summary.candidate==candidate].iloc[0];ci=intervals[candidate]
        guards=[]
        for family in ['K562','CD4','H1']:
            selected=screen[(screen.regime=='lineage_holdout')&(screen.family==family)]
            a=selected[selected.variant==candidate].iloc[0];b=selected[selected.variant=='identity_rank'].iloc[0]
            guards.append(bool(a.mse<=b.mse+1e-12 and a.fidelity_proxy>=b.fidelity_proxy-1e-12 and a.jaccard_proxy>=b.jaccard_proxy-1e-12))
        decisions[candidate]={'summary_context_guards':dict(zip(['K562','CD4','H1'],guards)),
            'fidelity_improvement_supported':ci['fidelity']['ci_low']>0,
            'overlap_improvement_supported':ci['jaccard']['ci_low']>0,
            'expression_improvement_supported':bool(model.expression_ratio<base.expression_ratio and ci['uncorrected_profile_mse']['ci_high']<0),
            'discrimination_preserved':bool(model.pds>=base.pds),
            'variance_preserved':bool(model.variance_ratio_to_template>=.9*base.variance_ratio_to_template)}
        decisions[candidate]['accepted']=all(guards) and all(value for key,value in decisions[candidate].items() if key not in ['summary_context_guards','accepted'])
    accepted=[c for c,d in decisions.items() if d['accepted'] and c!='identity_rank']
    selected=max(accepted,key=lambda c:float(summary.loc[summary.candidate==c,'fidelity'].iloc[0])) if accepted else 'v4'
    build.write_json(OUT/'adoption_decision.json',{'selected':selected,'decisions':decisions,'paired_target_bootstrap':intervals,
        'primary_targets':primary.tolist(),'all_source_targets':h1_names,'blend':.25,'seeds':2,'partitions':2,
        'source_count_sha256':extraction['counts_sha256'],'v4_components_sha256':record['source_components_sha256'],
        'count_script_sha256':build.sha256(Path(__file__)),'vcc_methods_sha256':build.sha256(ROOT/'v6_methods.py'),
        'head_methods_sha256':build.sha256(ROOT/'embedding_response_methods.py'),
        'baseline_scope':'Fixed V4 CD4/K562 transfer coefficients and C CD4 condition mix; all H1 direct/prior components removed. This is not the exact published A/B/C V4 artifact.',
        'evaluation_scope':'Development comparison on public H1 controls and cells, with H1 lineage excluded from fitting every head. Only 9 larger target groups expected. Metrics are unanchored, not a hidden leaderboard estimate.',
        'package_generated':False})
    build.write_json(OUT/'generated_cell_identity.json',identity)
    arrays={}
    for candidate in CANDIDATES:
        for part in [0,1]:
            for seed in [0,1]:
                for k,key in enumerate(['pred','truth','control','pred_jack','truth_jack','control_jack']):
                    arrays[f'{candidate}_{part}_{seed}_{key}']=np.stack([profiles[(str(t),part,seed,candidate)][k] for t in primary])
    np.savez_compressed(OUT/'generated_profiles.npz',genes=genes,targets=primary,**arrays)
    print(summary.to_string(index=False),flush=True);print('Adoption:',selected,flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
