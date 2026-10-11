"""V15 nested response, disjoint NTC, and generated-count validation."""
from pathlib import Path
import argparse,json,gc,hashlib
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
from v15_ensemble_methods import ROOT,A,FAMILIES,EXPERTS,cache_experts,fit_ensemble,predict
from v15_population_methods import source_quantiles,distribution_model,emit
from fit_v7_response import read_cache,raw_inputs
from v7_response_methods import summary_score
from v6_methods import moments
from v12_cpu_metrics import de,direction_metrics
from v8_methods import panel_metrics
from audit_transfer_validation import profile_jackknife_chunked,pds
from evaluate_atlas_transfer import load_atlas_model
import build_submission as build
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio']

def clean(v):
    if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [clean(x) for x in v]
    if isinstance(v,np.ndarray):return clean(v.tolist())
    if isinstance(v,(np.integer,)):return int(v)
    if isinstance(v,(float,np.floating)):return float(v) if np.isfinite(v) else None
    if isinstance(v,np.bool_):return bool(v)
    return v

def write(path,value):build.write_json(path,clean(value))
def load_plan():
    p=json.loads((A/'plan.json').read_text())
    for n,h in p['dependencies'].items():assert build.sha256(ROOT/n)==h,'Frozen V15 dependency changed: '+n
    old=json.loads((ROOT/'experiments/v13_response_programs/analysis/reference_model.json').read_text())
    for n,h in old['dependencies'].items():assert build.sha256(ROOT/n)==h,'V13 dependency changed: '+n
    assert build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')==p['source_cache_sha256']
    return p

def model(data,family):
    excluded=[] if family=='all' else [family];folder=A/'models';folder.mkdir(exist_ok=True);side=folder/(family+'.json')
    inputs={'plan':build.sha256(A/'plan.json'),'methods':build.sha256(ROOT/'v15_ensemble_methods.py')}
    if side.exists():
        record=json.loads(side.read_text());assert record['inputs']==inputs
        return cache_experts(data,excluded),record['model']
    base,record=fit_ensemble(data,excluded)
    if family!='all':
        frozen=ROOT/'experiments/v14_context_interactions/analysis/models'/(family+'.npz')
        with np.load(frozen,allow_pickle=False) as z:np.testing.assert_array_equal(base['v13'],z['base'])
    write(side,{'inputs':inputs,'model':record,'exact_v14_v13_bank_replay':family!='all'})
    return base,record

def screen():
    plan=load_plan();data=read_cache();rows=[];diversity=[]
    for family in FAMILIES:
        base,record=model(data,str(family));ids=np.flatnonzero((data['row_families']==family)&~np.isin(data['targets'][data['row_targets']],plan['reserved_targets']))
        for ci in np.unique(data['row_contexts'][ids]):
            bank,_=predict(data,base,record,data['profiles'][ci]);bank['target_shuffle']=bank['ensemble'][np.random.default_rng(2026101521).permutation(300)]
            for ri in ids[data['row_contexts'][ids]==ci]:
                ti=data['row_targets'][ri];mask=data['valid'][ri].copy();mask[base['own'][ti]]=False;truth=data['truth'][ri]
                for name,value in bank.items():rows.append({'family':str(family),'context':int(ci),'target':str(data['targets'][ti]),'candidate':name,**summary_score(value[ti],truth,data['errors'][ri],data['q'][ri],mask)})
                real=mask&(data['q'][ri]<.05)
                for name in ['student','specialist','ensemble']:
                    columns=np.flatnonzero(mask)[::10];a=bank['v13'][ti,columns]-truth[columns];b=bank[name][ti,columns]-truth[columns]
                    diversity.append({'family':str(family),'target':str(data['targets'][ti]),'expert':name,'error_correlation':float(np.corrcoef(a,b)[0,1]) if len(columns)>2 and a.std()>0 and b.std()>0 else np.nan,
                        'reference_DE_genes':int(real.sum()),'corrects_v13_sign':int((real&(np.sign(bank['v13'][ti])!=np.sign(truth))&(np.sign(bank[name][ti])==np.sign(truth))).sum()),
                        'breaks_v13_sign':int((real&(np.sign(bank['v13'][ti])==np.sign(truth))&(np.sign(bank[name][ti])!=np.sign(truth))).sum())})
        pd.DataFrame(rows).to_csv(A/'response_scores.csv',index=False);print('V15 outer response screen complete:',family,flush=True);del base;gc.collect()
    frame=pd.DataFrame(rows);fields=['mse_ratio','fidelity_proxy','jaccard_proxy','sign_reference_DE','selection_objective']
    family=frame.groupby(['candidate','family','context'])[fields].mean().groupby(['candidate','family']).mean();family.to_csv(A/'response_family.csv');summary=family.groupby('candidate').mean();summary.to_csv(A/'response_summary.csv')
    pd.DataFrame(diversity).to_csv(A/'expert_diversity.csv',index=False)
    winner=summary.loc[['student50','ensemble'],'selection_objective'].idxmax()
    result={'complete':True,'response':winner,'scores_sha256':build.sha256(A/'response_scores.csv'),'plan_sha256':build.sha256(A/'plan.json'),'scope':'Nested family-held calibration, previously viewed public source proxies. Counts and reserved algorithm verification required; no official score prediction.'}
    path=A/'response_selection.json'
    if path.exists():assert json.loads(path.read_text())==result
    else:write(path,result)
    print(summary.to_string(),flush=True)

def digest(x):
    h=hashlib.sha256()
    for a in [x.data,x.indices,x.indptr]:h.update(a.tobytes())
    return h.hexdigest()

def evaluate(stage,workers):
    plan=load_plan();data=read_cache();genes=data['genes'];targets=data['targets'];atlas,_=load_atlas_model()
    response=json.loads((A/'response_selection.json').read_text())['response']
    if stage=='null':names=['compact','quantile10','quantile25']
    elif stage=='development':names=plan['count_candidates']
    else:
        selection=json.loads((A/'selection.json').read_text());assert selection['candidate']!='v13'
        assert selection['scores_sha256']==build.sha256(A/'development_panels.csv');names=['v13',selection['candidate']]
    cache=A/(stage+'_cache');cache.mkdir(exist_ok=True);records=[];panels=[];cases=[]
    for ctx in plan['contexts']:
        family=ctx['family'];base,ensemble=model(data,family);x,meta,measured,_=raw_inputs(family,genes)
        population=x[ctx['left']];depth,p,b,_=moments(population);query=p*1e6;query[~measured]=np.nan
        bank,_=predict(data,base,ensemble,query);delta,supported,support=distribution_model(data,[family])
        ref=x[ctx['right']];cp,cj=profile_jackknife_chunked(ref)
        groups=ctx['null_cases'] if stage=='null' else [c for c in ctx['cells'] if c['split']==stage]
        byseed={}
        for cell in groups:byseed.setdefault(cell['seed'],[]).append(cell)
        for seed,cells in byseed.items():
            def one(cell):
                null=stage=='null';target='Non-Targeting' if null else cell['target'];ti=-1 if null else int(np.flatnonzero(targets==target)[0]);own=None if null else int(base['own'][ti])
                file=cache/f'{family}_{ti:03d}_{seed}.npz';side=file.with_suffix('.json')
                inputs={'plan':build.sha256(A/'plan.json'),'model':build.sha256(A/'models'/(family+'.json')),'response':response,'names':names}
                if file.exists() and side.exists():
                    rec=json.loads(side.read_text());assert rec['inputs']==inputs and rec['sha256']==build.sha256(file)
                    with np.load(file,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
                    return rec['records'],arrays,rec['provenance']
                actual=x[cell['actual']];template=x[cell['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64)
                assert len(cell['actual'])==len(set(cell['actual']))==len(libs);assert not set(cell['actual'])&set(cell['template'])
                tp,tj=profile_jackknife_chunked(actual);truth=de(actual,ref,genes);arrays={'truth':tp,'truth_jack':tj};local=[]
                for name in names:
                    mode=name if null else ('quantile10' if name.endswith('10') else 'quantile25' if name.endswith('25') else 'compact')
                    effect=np.zeros(len(genes)) if null else bank['v13' if name.startswith('v13') else response][ti]
                    d=np.zeros((2,len(genes),5),np.float32) if null else delta[ti];mask=measured&(p*1e6>5) if null else supported[ti]
                    if name=='shuffled_target25':
                        shuffled=np.random.default_rng(2026101523).permutation(300);d=delta[shuffled[ti]];mask=supported[shuffled[ti]]
                    count,diag=emit(template,np.arange(len(libs)),effect,own,seed+max(ti,0),mode,d,mask,prepared=(libs,p,b),atlas=atlas)
                    pp,pj=profile_jackknife_chunked(count);arrays[name+'_pred']=pp;arrays[name+'_jack']=pj
                    metrics=direction_metrics(de(count,ref,genes),truth,genes,-1 if null else own)
                    norm=count.multiply((1/np.asarray(count.sum(1)).ravel())[:,None]);variance=np.asarray(norm.power(2).mean(0)).ravel()-np.asarray(norm.mean(0)).ravel()**2
                    local.append({'family':family,'context':ctx['context'],'target':target,'seed':seed,'candidate':name,'reference_mode':'full_independent','actual_cells':actual.shape[0],'prediction_cells':count.shape[0],
                        'reference_cells':ref.shape[0],'count_sha256':digest(count),'total_CPM_variance':float(variance.sum()*1e12),**metrics,**diag})
                prov={'family':family,'target':target,'seed':seed,'actual_unique_cells':actual.shape[0],'template_unique_cells':len(libs),'response_training_excludes':family,'distribution_training_excludes':family}
                np.savez_compressed(file,**arrays);write(side,{'inputs':inputs,'sha256':build.sha256(file),'records':local,'provenance':prov})
                print('V15',stage,family,target,seed,'completed',flush=True);return clean(local),arrays,prov
            with ThreadPoolExecutor(max_workers=workers) as pool:outputs=list(pool.map(one,cells))
            for rr,_,prov in outputs:records.extend(rr);cases.append(prov)
            tp=np.asarray([a['truth'] for _,a,_ in outputs]);tj=np.asarray([a['truth_jack'] for _,a,_ in outputs])
            if stage!='null':
                mask=~np.isin(genes,targets);constant=np.broadcast_to(cp,tp.shape);np.testing.assert_allclose(pds(constant[:,mask],tp[:,mask],cp[mask]),.5,atol=1e-12)
            for name in names:
                pp=np.asarray([a[name+'_pred'] for _,a,_ in outputs]);pj=np.asarray([a[name+'_jack'] for _,a,_ in outputs])
                if stage=='null':
                    num=float(((pp-tp)**2-pj-tj).sum(1).mean());score={'common_expression_numerator':num}
                else:score=panel_metrics(pp,tp,cp,pj,tj,cj,genes,[c['target'] for c in cells],targets)
                subset=pd.DataFrame([r for rr,_,_ in outputs for r in rr if r['candidate']==name])
                panels.append({**subset.mean(numeric_only=True).to_dict(),'family':family,'context':ctx['context'],'seed':seed,'candidate':name,'reference_mode':'full_independent',**score})
            pd.DataFrame(records).to_csv(A/f'{stage}_scores.csv',index=False);pd.DataFrame(panels).to_csv(A/f'{stage}_panels.csv',index=False)
        del x,base,bank,delta;gc.collect()
    frame=pd.DataFrame(panels);fields=FIELDS if stage!='null' else ['predicted_calls','reference_calls','total_CPM_variance','common_expression_numerator']
    summary=frame.groupby(['candidate','family'])[fields].mean().groupby('candidate').mean();summary.to_csv(A/f'{stage}_summary.csv')
    write(A/f'{stage}_evaluation.json',{'complete':True,'population_cases':len(cases)*len(names),'target_seed_cases':len(cases),'provenance':cases,'scores_sha256':build.sha256(A/f'{stage}_scores.csv'),'panels_sha256':build.sha256(A/f'{stage}_panels.csv'),'plan_sha256':build.sha256(A/'plan.json')})
    print(summary.to_string(),flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['sources','screen','null','development','confirmation']);parser.add_argument('--workers',type=int,default=4);args=parser.parse_args()
    if not 1<=args.workers<=6:parser.error('workers must be 1..6')
    with threadpool_limits(limits=2):
        if args.stage=='sources':load_plan();source_quantiles(read_cache())
        elif args.stage=='screen':screen()
        else:evaluate(args.stage,args.workers)
