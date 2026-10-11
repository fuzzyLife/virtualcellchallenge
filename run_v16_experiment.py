"""V16 nested calibration, control-only diagnostics and paired count experiments."""
from pathlib import Path
import argparse,json,gc,hashlib
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
import build_submission as b
from v16_methods import ROOT,OUT,A,FAMILIES,CANDIDATES,source_cache,fit_calibration,predict,emit
from fit_v7_response import read_cache,raw_inputs
from v6_methods import moments
from v7_response_methods import summary_score
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked,pds
from v8_methods import panel_metrics
from evaluate_atlas_transfer import load_atlas_model
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio']

def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,np.ndarray):return clean(x.tolist())
    if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
    if isinstance(x,np.integer):return int(x)
    if isinstance(x,np.bool_):return bool(x)
    return x

def write(p,x):b.write_json(p,clean(x))
def digest(x):
    h=hashlib.sha256()
    for a in [x.data,x.indices,x.indptr]:h.update(a.tobytes())
    return h.hexdigest()

def load_plan():
    p=json.loads((A/'plan.json').read_text())
    for n,h in p['dependencies'].items():assert b.sha256(ROOT/n)==h,'Frozen V16 method changed: '+n
    assert b.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')==p['source_cache_sha256']
    assert b.sha256(ROOT/'experiments/v15_weighted_transport/analysis/plan.json')==p['parent_plan_sha256']
    old=json.loads((ROOT/'experiments/v13_response_programs/analysis/reference_model.json').read_text())
    for n,h in old['dependencies'].items():assert b.sha256(ROOT/n)==h,'V13 dependency changed: '+n
    return p

def model(data,family):
    plan=json.loads((A/'plan.json').read_text());excluded=[] if family=='all' else [family];bank=source_cache(data,excluded);folder=A/'models';folder.mkdir(exist_ok=True);file=folder/(family+'.json')
    inputs={'plan':b.sha256(A/'plan.json'),'methods':b.sha256(ROOT/'v16_methods.py')}
    if file.exists():
        r=json.loads(file.read_text());assert r['inputs']==inputs;return bank,r['model']
    m=fit_calibration(data,excluded,plan['reserved_targets']);write(file,{'inputs':inputs,'model':m});print('V16 calibration fitted:',family,flush=True);return bank,m

def fit_models():
    plan=load_plan();data=read_cache();records=[]
    for family in FAMILIES:
        bank,m=model(data,str(family));rows=np.flatnonzero((data['row_families']==family)&~np.isin(data['targets'][data['row_targets']],plan['reserved_targets']))
        for ci in np.unique(data['row_contexts'][rows]):
            pred,_=predict(bank,m,data['profiles'][ci])
            for ri in rows[data['row_contexts'][rows]==ci]:
                ti=data['row_targets'][ri];valid=data['valid'][ri].copy();valid[bank['own'][ti]]=False
                for name in ['v13','calibrated']:records.append({'family':str(family),'context':int(ci),'target':str(data['targets'][ti]),'candidate':name,**summary_score(pred[name][ti],data['truth'][ri],data['errors'][ri],data['q'][ri],valid)})
        del bank;gc.collect()
    frame=pd.DataFrame(records);frame.to_csv(A/'response_scores.csv',index=False)
    fields=['mse_ratio','fidelity_proxy','jaccard_proxy','selection_objective'];summary=frame.groupby(['candidate','family','context'])[fields].mean().groupby(['candidate','family']).mean().groupby('candidate').mean();summary.to_csv(A/'response_summary.csv')
    write(A/'models_complete.json',{'complete':True,'plan_sha256':b.sha256(A/'plan.json'),'scores_sha256':b.sha256(A/'response_scores.csv'),'scope':'Fixed calibration, no response selection from these proxies.'});print(summary.to_string(),flush=True)

def prior_case(ctx,cell,ti):
    parent=ROOT/'experiments/v15_weighted_transport/analysis';plan=json.loads((parent/'plan.json').read_text());c=next(v for v in plan['contexts'] if v['family']==ctx['family']);old=next(v for v in c['cells'] if v['target']==cell['target'] and v['seed']==cell['seed'] and v['split']==cell['split'])
    assert all(old[k]==cell[k] for k in ['actual','template','reference']) and c['left']==ctx['left'] and c['right']==ctx['right']
    replay=json.loads((parent/'same_input_baseline_replay.json').read_text());assert replay['passed'] and replay['cases']==18
    file=parent/'development_cache'/f'{ctx["family"]}_{ti:03d}_{cell["seed"]}.npz';rec=json.loads(file.with_suffix('.json').read_text());assert rec['sha256']==b.sha256(file)
    with np.load(file,allow_pickle=False) as z:arrays={k:z[k].copy() for k in ['truth','truth_jack','v13_pred','v13_jack']}
    record=next(r.copy() for r in rec['records'] if r['candidate']=='v13');record['context']=ctx['context'];record['reused_v15_baseline']=True
    return record,arrays,rec['sha256']

def evaluate(stage,workers):
    plan=load_plan();data=read_cache();genes=data['genes'];targets=data['targets'];atlas,_=load_atlas_model()
    selection=json.loads((A/'selection.json').read_text()) if stage=='confirmation' else None
    if stage=='null':names=['compact','balanced','residual10']
    elif stage=='development':names=list(CANDIDATES)
    else:
        assert selection['candidate']!='v13';names=['v13',selection['candidate']]
        if selection['candidate'].startswith('calibrated_'):names.append(selection['candidate'].replace('calibrated_','uniform_'))
    cache=A/(stage+'_cache');cache.mkdir(exist_ok=True);scores=[];panels=[];provenance=[]
    for ctx in plan['contexts']:
        family=ctx['family'];x,meta,measured,_=raw_inputs(family,genes);population=x[ctx['left']];_,p,bulk,_=moments(population);query=p*1e6;query[~measured]=np.nan
        if stage!='null':bank,m=model(data,family);responses,factors=predict(bank,m,query);model_hash=b.sha256(A/'models'/(family+'.json'))
        else:responses=None;model_hash=None
        ref=x[ctx['right']];cp,cj=profile_jackknife_chunked(ref);_,rp,_,rv=moments(ref);rd=np.asarray((ref>0).mean(0)).ravel();eligible=rp*1e6>5
        cases=ctx['null_cases'] if stage=='null' else [c for c in ctx['cells'] if c['split']==stage]
        byseed={}
        for c in cases:byseed.setdefault(c['seed'],[]).append(c)
        for seed,cells in byseed.items():
            def one(cell):
                null=stage=='null';target='Non-Targeting' if null else cell['target'];ti=-1 if null else int(np.flatnonzero(targets==target)[0]);own=None if null else int(bank['own'][ti]);file=cache/f'{family}_{ti:03d}_{seed}.npz';side=file.with_suffix('.json')
                inputs={'plan':b.sha256(A/'plan.json'),'model':model_hash,'names':names}
                if side.exists():
                    rec=json.loads(side.read_text());assert rec['inputs']==inputs and rec['sha256']==b.sha256(file)
                    with np.load(file,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
                    return rec['records'],arrays,rec['provenance']
                actual=x[cell['actual']];template=x[cell['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64);assert len(set(cell['template']))==400 and not set(cell['actual'])&set(cell['template'])
                truth=de(actual,ref,genes);tp,tj=profile_jackknife_chunked(actual);arrays={'truth':tp,'truth_jack':tj};local=[];reused=0
                for name in names:
                    if name=='v13' and stage=='development':
                        r,a,prior_sha=prior_case(ctx,cell,ti);np.testing.assert_array_equal(tp,a['truth']);np.testing.assert_array_equal(tj,a['truth_jack']);arrays.update(a);local.append(r);reused+=1;continue
                    if null:kind=name;effect=np.zeros(len(genes))
                    else:response,kind=CANDIDATES[name];effect=responses[response][ti]
                    count,diag=emit(template,np.arange(400),effect,own,seed+max(ti,0),kind,prepared=(libs,p,bulk),atlas=atlas)
                    pp,pj=profile_jackknife_chunked(count);arrays[name+'_pred']=pp;arrays[name+'_jack']=pj;metrics=direction_metrics(de(count,ref,genes),truth,genes,-1 if own is None else own)
                    norm=count.multiply((1/libs)[:,None]);var=np.asarray(norm.power(2).mean(0)).ravel()-np.asarray(norm.mean(0)).ravel()**2;detect=np.asarray((count>0).mean(0)).ravel()
                    local.append({'family':family,'context':ctx['context'],'target':target,'seed':seed,'candidate':name,'actual_cells':actual.shape[0],'prediction_cells':400,'reference_cells':ref.shape[0],'count_sha256':digest(count),'reference_mode':'full_independent',
                        'reused_v15_baseline':False,'variance_ratio':float(var.sum()/rv.sum()),'variance_log_distance':float(abs(np.log(max(var.sum()/rv.sum(),1e-20)))),'detection_MAE':float(abs(detect-rd)[eligible].mean()),**metrics,**diag})
                prov={'family':family,'target':target,'seed':seed,'actual_unique_cells':actual.shape[0],'template_unique_cells':400,'held_family':family,'compared_populations':len(names),'newly_generated_populations':len(names)-reused,'reused_populations':reused}
                np.savez_compressed(file,**arrays);write(side,{'inputs':inputs,'sha256':b.sha256(file),'records':local,'provenance':prov});print('V16',stage,family,target,seed,'completed',flush=True);return clean(local),arrays,prov
            with ThreadPoolExecutor(max_workers=workers) as pool:outputs=list(pool.map(one,cells))
            for rr,_,prov in outputs:scores.extend(rr);provenance.append(prov)
            tp=np.asarray([a['truth'] for _,a,_ in outputs]);tj=np.asarray([a['truth_jack'] for _,a,_ in outputs])
            if stage!='null':
                axis=~np.isin(genes,targets);np.testing.assert_allclose(pds(np.broadcast_to(cp,tp.shape)[:,axis],tp[:,axis],cp[axis]),.5,atol=1e-12)
            for name in names:
                pp=np.asarray([a[name+'_pred'] for _,a,_ in outputs]);pj=np.asarray([a[name+'_jack'] for _,a,_ in outputs]);metric={} if stage=='null' else panel_metrics(pp,tp,cp,pj,tj,cj,genes,[c['target'] for c in cells],targets)
                f=pd.DataFrame([r for rr,_,_ in outputs for r in rr if r['candidate']==name]);panels.append({**f.mean(numeric_only=True).to_dict(),'family':family,'context':ctx['context'],'seed':seed,'candidate':name,**metric})
            pd.DataFrame(scores).to_csv(A/f'{stage}_scores.csv',index=False);pd.DataFrame(panels).to_csv(A/f'{stage}_panels.csv',index=False)
        del x,population;gc.collect()
    frame=pd.DataFrame(panels);fields=['predicted_calls','reference_calls','variance_ratio','variance_log_distance','detection_MAE'] if stage=='null' else FIELDS
    summary=frame.groupby(['candidate','family'])[fields].mean().groupby('candidate').mean();summary.to_csv(A/f'{stage}_summary.csv')
    write(A/f'{stage}_evaluation.json',{'complete':True,'population_comparisons':len(scores),'newly_generated_populations':sum(v['newly_generated_populations'] for v in provenance),'reused_populations':sum(v['reused_populations'] for v in provenance),'target_seed_cases':len(provenance),'provenance':provenance,'scores_sha256':b.sha256(A/f'{stage}_scores.csv'),'panels_sha256':b.sha256(A/f'{stage}_panels.csv'),'plan_sha256':b.sha256(A/'plan.json')});print(summary.to_string(),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['models','null','development','confirmation']);p.add_argument('--workers',type=int,default=6);args=p.parse_args()
    if not 1<=args.workers<=6:p.error('workers must be 1..6')
    with threadpool_limits(limits=2):
        if args.stage=='models':fit_models()
        else:evaluate(args.stage,args.workers)
