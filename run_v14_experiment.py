"""V14 response screen and independent-treatment oracle/count experiments."""
from pathlib import Path
import argparse,gc,json,hashlib
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
from v14_context_methods import A,ROOT,NAMES,FAMILIES,fit,predict,family_profiles
from fit_v7_response import read_cache,raw_inputs
from v7_response_methods import summary_score
from v6_methods import moments
from v11_methods import emit,SETTINGS
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked,pds
from v8_methods import panel_metrics
from evaluate_atlas_transfer import load_atlas_model
import build_submission as build
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio']


def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,(np.integer,)):return int(x)
    if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
    if isinstance(x,np.bool_):return bool(x)
    return x


def write(path,data):build.write_json(path,clean(data))


def load_plan():
    plan=json.loads((A/'plan.json').read_text())
    for name,digest in plan['dependencies'].items():assert build.sha256(ROOT/name)==digest,'Frozen dependency changed: '+name
    assert build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')==plan['source_cache_sha256']
    freeze=A/'execution_freeze.json'
    files=['run_v14_experiment.py','test_v14_context.py']
    record={'plan_sha256':build.sha256(A/'plan.json'),'dependencies':{f:build.sha256(ROOT/f) for f in files}}
    if freeze.exists():assert json.loads(freeze.read_text())==record,'Execution code changed after freeze'
    else:write(freeze,record)
    return plan


def model_cache(data,family):
    folder=A/'models';folder.mkdir(exist_ok=True);file=folder/(family+'.npz');side=file.with_suffix('.json')
    inputs={'plan':build.sha256(A/'plan.json'),'source':build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),'methods':build.sha256(ROOT/'v14_context_methods.py')}
    if file.exists() and side.exists():
        record=json.loads(side.read_text());assert record['inputs']==inputs and record['sha256']==build.sha256(file)
        with np.load(file,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
        return {'base':arrays['base'],'own':arrays['own'],'features':{k:arrays['f_'+k] for k in ['axis','center','basis','scale','radius']},
                'slopes':{k:arrays['s_'+k] for k in ['context2','context8','shuffled2']},'centers':{k:arrays['c_'+k] for k in ['context2','context8','shuffled2']},**record['metadata']}
    model=fit(data,family)
    arrays={'base':model['base'],'own':model['own'],**{'f_'+k:v for k,v in model['features'].items()},
            **{'s_'+k:v for k,v in model['slopes'].items()},**{'c_'+k:v for k,v in model['centers'].items()}}
    np.savez_compressed(file,**arrays)
    write(side,{'inputs':inputs,'sha256':build.sha256(file),'metadata':{k:model[k] for k in ['excluded','training_families','training_rows','feature_genes']}})
    return model


def balanced(frame):
    return frame.groupby(['candidate','family'])[FIELDS].mean().groupby('candidate').mean().reset_index()


def screen():
    plan=load_plan();data=read_cache();records=[];features=[]
    for family in FAMILIES:
        model=model_cache(data,str(family))
        rows=np.flatnonzero((data['row_families']==family)&~np.isin(data['targets'][data['row_targets']],plan['confirmation_targets']))
        for ci in np.unique(data['row_contexts'][rows]):
            bank,diag=predict(model,data['profiles'][ci]);features.append({'family':str(family),'context':int(ci),**diag,'feature_genes':model['feature_genes']})
            bank['zero']=np.zeros_like(bank['v13'])
            bank['generic']=np.broadcast_to(bank['v13'].mean(0),bank['v13'].shape)
            bank['target_shuffle']=bank['v13'][np.random.default_rng(2026101417).permutation(len(data['targets']))]
            for ri in rows[data['row_contexts'][rows]==ci]:
                ti=data['row_targets'][ri];mask=data['valid'][ri].copy();mask[model['own'][ti]]=False
                for name,value in bank.items():
                    metrics=summary_score(value[ti],data['truth'][ri],data['errors'][ri],data['q'][ri],mask)
                    records.append({'family':str(family),'context':int(ci),'target':str(data['targets'][ti]),'candidate':name,**metrics})
        pd.DataFrame(records).to_csv(A/'screen_scores.csv',index=False)
        print('V14 response screen completed',family,flush=True);del model;gc.collect()
    fields=['mse_ratio','fidelity_proxy','jaccard_proxy','sign_reference_DE','selection_objective']
    frame=pd.DataFrame(records);family=frame.groupby(['candidate','family','context'])[fields].mean().groupby(['candidate','family']).mean()
    family.to_csv(A/'screen_family.csv');summary=family.groupby('candidate').mean();summary.to_csv(A/'screen_summary.csv')
    winner=summary.loc[['context2','context8'],'selection_objective'].idxmax()
    write(A/'screen_result.json',{'completed':True,'advance':winner,'count_candidates':['v13',winner,'shuffled2','amplitude110'],
          'feature_diagnostics':features,'plan_sha256':build.sha256(A/'plan.json'),'scores_sha256':build.sha256(A/'screen_scores.csv'),
          'scope':'Previously viewed public summary proxies; count validation required. Whole evaluated family excluded; control feature basis trained only on other families. Reserved count-confirmation targets omitted from summary evaluation.'})
    print(summary.to_string(),flush=True)


def count_digest(x):
    h=hashlib.sha256()
    for a in [x.data,x.indices,x.indptr]:h.update(a.tobytes())
    return h.hexdigest()


def compact_oracle(libraries,p,b,seed,atlas):
    rng=np.random.default_rng(seed);expected=libraries[:,None]*p[None,:]
    template=.35*expected+.65*rng.poisson(expected);template/=template.sum(1,keepdims=True)
    counts=atlas.dual_moment_counts(template,p,b,depths=libraries,seed=seed,iterations=200)
    np.testing.assert_array_equal(counts.sum(1,dtype=np.int64),libraries)
    return sparse.csr_matrix(counts)


def evaluate(stage,workers):
    plan=load_plan();data=read_cache();genes=data['genes'];targets=data['targets'];atlas,_=load_atlas_model()
    if stage=='oracle':names=['v13_compact','oracle_moments_compact','actual_treatment_replica']
    elif stage=='development':names=json.loads((A/'screen_result.json').read_text())['count_candidates']
    else:
        frozen=json.loads((A/'selection.json').read_text());assert frozen['promote'] and not (A/'confirmation_evaluation.json').exists()
        assert frozen['development_sha256']==build.sha256(A/'development_summary.csv')
        names=['v13',frozen['candidate']]
    cache=A/(stage+'_cache');cache.mkdir(exist_ok=True);records=[];panels=[];provenance=[]
    for ctx in plan['contexts']:
        family=ctx['family'];model=model_cache(data,family);x,meta,measured,_=raw_inputs(family,genes)
        population=x[ctx['left']];depth,p,b,var=moments(population);query=p*1e6;query[~measured]=np.nan
        bank,feature_diag=predict(model,query)
        assert np.isfinite(np.array(list(bank.values()))).all()
        refs={'full_independent':x[ctx['right']]}
        if stage=='oracle':
            groups=[c for c in plan['oracle_cases'] if c['family']==family]
        else:groups=[c for c in ctx['cells'] if c['split']==stage]
        if not groups:del x,model;gc.collect();continue
        byseed={}
        for cell in groups:byseed.setdefault(cell['seed'],[]).append(cell)
        for seed,cells in byseed.items():
            if stage!='oracle':refs['sampled400']=x[cells[0]['reference']]
            reference_profiles={mode:profile_jackknife_chunked(ref) for mode,ref in refs.items()}
            def one(cell):
                target=cell['target'];ti=int(np.flatnonzero(targets==target)[0]);own=int(model['own'][ti]);file=cache/f'{family}_{ti:03d}_{seed}.npz';side=file.with_suffix('.json')
                inputs={'plan':build.sha256(A/'plan.json'),'execution':build.sha256(A/'execution_freeze.json'),'model':build.sha256(A/'models'/f'{family}.npz'),'names':names}
                if file.exists() and side.exists():
                    record=json.loads(side.read_text());assert record['inputs']==inputs and record['sha256']==build.sha256(file)
                    with np.load(file,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
                    return record['records'],arrays,record['provenance']
                if stage=='oracle':
                    calibration=x[cell['calibration']];actual=x[cell['truth']]
                    assert not set(cell['calibration'])&set(cell['truth']) and calibration.shape[0]==actual.shape[0]==200
                    libraries,op,ob,_=moments(calibration);template=calibration;libs=libraries.astype(np.int64)
                    provenance={'family':family,'target':target,'seed':seed,'oracle_calibration_cells':200,'truth_cells':200,
                                'shared_depths':'Known treated-calibration depths for all oracle arms, diagnostic only; distinct from deployable count tests.'}
                else:
                    actual=x[cell['actual']];template=x[cell['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64)
                    assert actual.shape[0]==400 and len(set(cell['actual']))==400
                    provenance={'family':family,'target':target,'seed':seed,'truth_cells':400,'oracle_calibration_cells':0,'shared_depths':'Actual disjoint template-control libraries.'}
                tp,tj=profile_jackknife_chunked(actual);truth={mode:de(actual,ref,genes) for mode,ref in refs.items()};arrays={'truth':tp,'truth_jack':tj};rows=[]
                for name in names:
                    if name=='actual_treatment_replica':count=calibration.copy();diag={}
                    elif name=='oracle_moments_compact':count=compact_oracle(libs,op,ob,seed+ti,atlas);diag={}
                    else:
                        effect=bank['v13' if name=='v13_compact' else name][ti]
                        count,diag=emit(template,np.arange(len(libs)),effect,own,seed+ti,SETTINGS['consensus20_noise65'],prepared=(libs,p,b),atlas=atlas)
                    np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),libs)
                    assert np.issubdtype(count.dtype,np.integer) and (count.data>=0).all()
                    pp,pj=profile_jackknife_chunked(count);arrays[name+'_pred']=pp;arrays[name+'_jack']=pj
                    for mode,ref in refs.items():
                        metrics=direction_metrics(de(count,ref,genes),truth[mode],genes,own)
                        rows.append({'family':family,'context':ctx['context'],'target':target,'seed':seed,'stage':stage,'candidate':name,'reference_mode':mode,
                                     'actual_cells':actual.shape[0],'prediction_cells':count.shape[0],'reference_cells':ref.shape[0],
                                     'count_sha256':count_digest(count),**metrics,**diag})
                    del count
                np.savez_compressed(file,**arrays);write(side,{'inputs':inputs,'sha256':build.sha256(file),'records':rows,'provenance':provenance})
                print('V14',stage,family,target,seed,'completed',flush=True)
                return clean(rows),arrays,provenance
            with ThreadPoolExecutor(max_workers=workers) as executor:outputs=list(executor.map(one,cells))
            for rows,_,prov in outputs:records.extend(rows);provenance.append(prov)
            tp=np.array([a['truth'] for _,a,_ in outputs]);tj=np.array([a['truth_jack'] for _,a,_ in outputs])
            # Oracle seeds vary by target. Profile panels are pooled below after all family cases.
            if stage!='oracle':
                mask=~np.isin(genes,targets);constant=np.broadcast_to(reference_profiles['full_independent'][0],tp.shape)
                np.testing.assert_allclose(pds(constant[:,mask],tp[:,mask],reference_profiles['full_independent'][0][mask]),.5,atol=1e-12)
                for mode,(cp,cj) in reference_profiles.items():
                    for name in names:
                        pp=np.array([a[name+'_pred'] for _,a,_ in outputs]);pj=np.array([a[name+'_jack'] for _,a,_ in outputs])
                        score=panel_metrics(pp,tp,cp,pj,tj,cj,genes,[c['target'] for c in cells],targets)
                        subset=pd.DataFrame([r for rows,_,_ in outputs for r in rows if r['candidate']==name and r['reference_mode']==mode])
                        panels.append({'family':family,'context':ctx['context'],'seed':seed,'candidate':name,'reference_mode':mode,**subset.mean(numeric_only=True).to_dict(),**score})
            pd.DataFrame(records).to_csv(A/f'{stage}_scores.csv',index=False)
        if stage=='oracle':
            combined=[]
            for cell in groups:
                ti=int(np.flatnonzero(targets==cell['target'])[0]);file=cache/f'{family}_{ti:03d}_{cell["seed"]}.npz'
                with np.load(file,allow_pickle=False) as z:combined.append({k:z[k].copy() for k in z.files})
            tp=np.array([a['truth'] for a in combined]);tj=np.array([a['truth_jack'] for a in combined]);cp,cj=reference_profiles['full_independent']
            for name in names:
                pp=np.array([a[name+'_pred'] for a in combined]);pj=np.array([a[name+'_jack'] for a in combined])
                score=panel_metrics(pp,tp,cp,pj,tj,cj,genes,[c['target'] for c in groups],targets)
                subset=pd.DataFrame([r for r in records if r['family']==family and r['candidate']==name])
                panels.append({'family':family,'context':ctx['context'],'seed':-1,'candidate':name,'reference_mode':'full_independent',**subset.mean(numeric_only=True).to_dict(),**score})
        pd.DataFrame(panels).to_csv(A/f'{stage}_panels.csv',index=False)
        summary=balanced(pd.DataFrame(panels)[pd.DataFrame(panels).reference_mode=='full_independent']);summary.to_csv(A/f'{stage}_summary.csv',index=False)
        del x,model,bank;gc.collect()
    write(A/f'{stage}_evaluation.json',{'completed':True,'names':names,'generated_or_real_population_cases':len(provenance)*len(names),'target_seed_cases':len(provenance),
          'provenance':provenance,'plan_sha256':build.sha256(A/'plan.json'),'scores_sha256':build.sha256(A/f'{stage}_scores.csv'),'summary_sha256':build.sha256(A/f'{stage}_summary.csv'),
          'scope':'Unanchored public counts; oracle arms use held treatments only for diagnostic ceilings. Development/confirmation responses exclude the entire scored family.'})
    print(pd.read_csv(A/f'{stage}_summary.csv').to_string(index=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['screen','oracle','development','confirmation'])
    parser.add_argument('--workers',type=int,default=3)
    args=parser.parse_args()
    if not 1<=args.workers<=6:parser.error('workers must be 1..6')
    with threadpool_limits(limits=2):
        if args.stage=='screen':screen()
        else:evaluate(args.stage,args.workers)
