"""Paired complete-count V21 public-family experiment with independent NTC references."""
import argparse,json,gc,hashlib
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
import build_submission as b
from v21_methods import ROOT,A,NAMES,FAMILIES,load_plan,load_query,load_expert,target_effects,emit
from fit_v7_response import read_cache,raw_inputs
from v16_methods import source_cache
from v8_methods import gated,panel_metrics
from v6_methods import moments
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked,pds
from evaluate_atlas_transfer import load_atlas_model
from scipy.stats import wasserstein_distance
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio','shape_wasserstein','covariance_error','detection_error']

def program_scores(count,query_inputs):
 depth=np.asarray(count.sum(1)).ravel();x=count[:,query_inputs['positions']].astype(float).multiply((1e4/depth)[:,None]).tocsr();x.data=np.log1p(x.data);z=np.asarray(x@(query_inputs['basis']/query_inputs['scale'][:,None]));return z-z.mean(0)
def clean(x):
 if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
 if isinstance(x,(list,tuple,np.ndarray)):return [clean(v) for v in x]
 if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
 if isinstance(x,np.integer):return int(x)
 if isinstance(x,np.bool_):return bool(x)
 return x
def write(p,x):b.write_json(p,clean(x))
def digest(x):
 h=hashlib.sha256()
 for a in [x.data,x.indices,x.indptr]:h.update(a.tobytes())
 return h.hexdigest()
def evaluate(stage,workers):
 plan=load_plan();data=read_cache();genes=data['genes'];targets=data['targets'];atlas,_=load_atlas_model();scores=[];panels=[];provenance=[];diagnostics=[]
 if stage=='confirmation':assert json.loads((A/'selection.json').read_text())['new_candidate_passes']
 names=NAMES if stage=='development' else ['v16',json.loads((A/'selection.json').read_text())['candidate'],'amplitude','shuffled_joint','stack_mean','reweight'];names=list(dict.fromkeys(names))
 cache=A/(stage+'_cache');cache.mkdir(exist_ok=True)
 for ctx in plan['contexts']:
  family=ctx['family'];x,meta,measured,_=raw_inputs(family,genes);population=x[ctx['left']];_,p,bulk,_=moments(population);query=p*1e6;query[~measured]=np.nan;ref=x[ctx['right']];cp,cj=profile_jackknife_chunked(ref)
  bank=source_cache(data,[family]);responses=np.asarray([gated(v.copy(),np.nan_to_num(query),int(o)) for v,o in zip(bank['base'],bank['own'])],np.float32);query_inputs=load_query(family);expert=load_expert(family,stage);assert np.array_equal(expert['genes'],genes) and np.array_equal(expert['targets'],targets)
  byseed={}
  for c in [r for r in ctx['cells'] if r['split']==stage]:byseed.setdefault(c['seed'],[]).append(c)
  for seed,cells in byseed.items():
   def one(cell):
    target=cell['target'];ti=int(np.flatnonzero(targets==target)[0]);own=int(bank['own'][ti]);effect=responses[ti];file=cache/f'{family}_{target}_{seed}.npz';side=file.with_suffix('.json');inputs={'plan':b.sha256(A/'plan.json'),'expert':b.sha256(A/'experts'/f'{family}_{stage}.npz'),'query':b.sha256(A/'inputs'/f'{family}_query.npz'),'names':names}
    if side.exists():
     rec=json.loads(side.read_text());assert rec['inputs']==inputs and rec['sha256']==b.sha256(file)
     with np.load(file,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
     return rec['records'],arrays,rec['provenance'],rec['diagnostic']
    assert expert['valid'][ti];effects=target_effects(effect,expert,ti,query,own)
    actual=x[cell['actual']];template=x[cell['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64);assert len(set(cell['template']))==400 and not set(cell['actual'])&set(cell['template']) and not set(cell['template'])&set(ctx['right'])
    truth=de(actual,ref,genes);tp,tj=profile_jackknife_chunked(actual);actual_shape=program_scores(actual,query_inputs);directions=np.random.default_rng(20020).normal(size=(16,32));directions/=np.linalg.norm(directions,axis=0);truth_proj=actual_shape@directions;truth_cov=np.cov(actual_shape,rowvar=False);detection=np.asarray((actual>0).mean(0)).ravel();shape_mask=measured&~np.isin(genes,targets)&(np.nan_to_num(query)>5);arrays={'truth':tp,'truth_jack':tj};local=[];calls={}
    for name in names:
     count,diag=emit(template,effects[name],own,seed+ti,name,prepared=(libs,p,bulk),atlas=atlas,expert=expert,ti=ti,query_inputs=query_inputs,template_rows=cell['template'])
     pp,pj=profile_jackknife_chunked(count);arrays[name+'_pred']=pp;arrays[name+'_jack']=pj;stat=de(count,ref,genes);metrics=direction_metrics(stat,truth,genes,own);calls[name]=stat;pred_shape=program_scores(count,query_inputs);pred_proj=pred_shape@directions;shape=float(np.mean([wasserstein_distance(truth_proj[:,j],pred_proj[:,j]) for j in range(32)]));covariance=float(np.linalg.norm(np.cov(pred_shape,rowvar=False)-truth_cov)/max(np.linalg.norm(truth_cov),1e-12));detect=float(abs(np.asarray((count>0).mean(0)).ravel()[shape_mask]-detection[shape_mask]).mean())
     local.append({'family':family,'context':ctx['context'],'target':target,'seed':seed,'candidate':name,'actual_cells':actual.shape[0],'prediction_cells':400,'reference_cells':ref.shape[0],'count_sha256':digest(count),'shape_wasserstein':shape,'covariance_error':covariance,'detection_error':detect,**metrics,**diag})
    mask=truth['eligible']&(np.arange(len(genes))!=own);error0=effect[mask]-truth['lfc'][mask];error1=effects['stack_mean'][mask]-truth['lfc'][mask];true=mask&(truth['q']<.05);basecorrect=(calls['v16']['q']<.05)&(np.sign(calls['v16']['lfc'])==np.sign(truth['lfc']));chosen='stack_mean';newcorrect=(calls[chosen]['q']<.05)&(np.sign(calls[chosen]['lfc'])==np.sign(truth['lfc']))
    diagnostic={'family':family,'target':target,'seed':seed,'eligible_genes':int(mask.sum()),'baseline_expected_LFC_error':float(abs(error0).mean()),'Stack_expected_LFC_error':float(abs(error1).mean()),'changed_coordinates':int(np.count_nonzero(effects['stack_mean']!=effect)),'corrects_true_DE':int((true&~basecorrect&newcorrect).sum()),'breaks_true_DE':int((true&basecorrect&~newcorrect).sum())}
    prov={'family':family,'target':target,'seed':seed,'held_family':family,'actual_unique_cells':actual.shape[0],'template_unique_cells':400,'compared_populations':len(names)};np.savez_compressed(file,**arrays);write(side,{'inputs':inputs,'sha256':b.sha256(file),'records':local,'provenance':prov,'diagnostic':diagnostic});print('V21',stage,family,target,seed,'completed',flush=True);return clean(local),arrays,prov,diagnostic
   with ThreadPoolExecutor(max_workers=workers) as pool:outputs=list(pool.map(one,cells))
   for rr,_,prov,diag in outputs:scores.extend(rr);provenance.append(prov);diagnostics.append(diag)
   tp=np.asarray([a['truth'] for _,a,_,_ in outputs]);tj=np.asarray([a['truth_jack'] for _,a,_,_ in outputs]);axis=~np.isin(genes,targets);np.testing.assert_allclose(pds(np.broadcast_to(cp,tp.shape)[:,axis],tp[:,axis],cp[axis]),.5,atol=1e-12)
   for name in names:
    pp=np.asarray([a[name+'_pred'] for _,a,_,_ in outputs]);pj=np.asarray([a[name+'_jack'] for _,a,_,_ in outputs]);metric=panel_metrics(pp,tp,cp,pj,tj,cj,genes,[c['target'] for c in cells],targets);f=pd.DataFrame([r for rr,_,_,_ in outputs for r in rr if r['candidate']==name]);panels.append({**f.mean(numeric_only=True).to_dict(),'family':family,'context':ctx['context'],'seed':seed,'candidate':name,**metric})
   pd.DataFrame(scores).to_csv(A/f'{stage}_scores.csv',index=False);pd.DataFrame(panels).to_csv(A/f'{stage}_panels.csv',index=False);pd.DataFrame(diagnostics).to_csv(A/f'{stage}_expert_diagnostics.csv',index=False)
  del x,population;gc.collect()
 frame=pd.DataFrame(panels);summary=frame.groupby(['candidate','family'])[FIELDS].mean().groupby('candidate').mean();summary.to_csv(A/f'{stage}_summary.csv');write(A/f'{stage}_evaluation.json',{'complete':True,'population_comparisons':len(scores),'target_seed_cases':len(provenance),'provenance':provenance,'scores_sha256':b.sha256(A/f'{stage}_scores.csv'),'panels_sha256':b.sha256(A/f'{stage}_panels.csv'),'plan_sha256':b.sha256(A/'plan.json')});print(summary.to_string(),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);p.add_argument('--workers',type=int,default=6);a=p.parse_args();assert 1<=a.workers<=6
 with threadpool_limits(limits=2):evaluate(a.stage,a.workers)
