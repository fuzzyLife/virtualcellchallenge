"""Paired complete-count V18 public-family experiment with independent NTC references."""
import argparse,json,gc,hashlib
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
import build_submission as b
from v18_methods import ROOT,A,NAMES,load_plan,load_predictions,blend
from fit_v7_response import read_cache,raw_inputs
from v16_methods import source_cache,emit
from v8_methods import gated,panel_metrics
from v6_methods import moments
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked,pds
from evaluate_atlas_transfer import load_atlas_model
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio']
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
 names=NAMES if stage=='development' else ['v16',json.loads((A/'selection.json').read_text())['candidate'],'amplitude','shuffled'];names=list(dict.fromkeys(names))
 cache=A/(stage+'_cache');cache.mkdir(exist_ok=True)
 for ctx in plan['contexts']:
  family=ctx['family'];x,meta,measured,_=raw_inputs(family,genes);population=x[ctx['left']];_,p,bulk,_=moments(population);query=p*1e6;query[~measured]=np.nan;ref=x[ctx['right']];cp,cj=profile_jackknife_chunked(ref)
  bank=source_cache(data,[family]);responses=np.asarray([gated(v.copy(),np.nan_to_num(query),int(o)) for v,o in zip(bank['base'],bank['own'])],np.float32);expert=load_predictions(family,stage);positions=expert['positions'];support=np.zeros(len(genes),bool);support[positions]=expert['measured'];assert np.array_equal(genes[positions],expert['genes']) and np.array_equal(targets,expert['targets'])
  byseed={}
  for c in [r for r in ctx['cells'] if r['split']==stage]:byseed.setdefault(c['seed'],[]).append(c)
  for seed,cells in byseed.items():
   def one(cell):
    target=cell['target'];ti=int(np.flatnonzero(targets==target)[0]);own=int(bank['own'][ti]);effect=responses[ti];file=cache/f'{family}_{target}_{seed}.npz';side=file.with_suffix('.json');inputs={'plan':b.sha256(A/'plan.json'),'expert':b.sha256(A/'inference'/f'{family}_{stage}.npz'),'names':names}
    if side.exists():
     rec=json.loads(side.read_text());assert rec['inputs']==inputs and rec['sha256']==b.sha256(file)
     with np.load(file,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
     return rec['records'],arrays,rec['provenance'],rec['diagnostic']
    def mapped(i):
     assert expert['valid'][i];lfc=np.zeros(len(genes),np.float32);prob=lfc.copy();lfc[positions]=np.clip(expert['lfc'][i],-3,3);prob[positions]=expert['p_de'][i];return gated(lfc,np.nan_to_num(query),own),prob
    lf,prob=mapped(ti);effects=blend(effect,lf,prob,support,own);sl,sp=mapped(plan['shuffle'][ti]);effects['shuffled']=blend(effect,sl,sp,support,own)['confidence']
    actual=x[cell['actual']];template=x[cell['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64);assert len(set(cell['template']))==400 and not set(cell['actual'])&set(cell['template']) and not set(cell['template'])&set(ctx['right'])
    truth=de(actual,ref,genes);tp,tj=profile_jackknife_chunked(actual);arrays={'truth':tp,'truth_jack':tj};local=[];calls={}
    for name in names:
     count,diag=emit(template,np.arange(400),effects[name],own,seed+ti,'residual10',prepared=(libs,p,bulk),atlas=atlas)
     pp,pj=profile_jackknife_chunked(count);arrays[name+'_pred']=pp;arrays[name+'_jack']=pj;stat=de(count,ref,genes);metrics=direction_metrics(stat,truth,genes,own);calls[name]=stat
     local.append({'family':family,'context':ctx['context'],'target':target,'seed':seed,'candidate':name,'actual_cells':actual.shape[0],'prediction_cells':400,'reference_cells':ref.shape[0],'count_sha256':digest(count),**metrics,**diag})
    mask=truth['eligible']&support&(np.arange(len(genes))!=own);error0=effect[mask]-truth['lfc'][mask];error1=lf[mask]-truth['lfc'][mask];cor=float(np.corrcoef(error0,error1)[0,1]) if mask.sum()>2 else None
    true=mask&(truth['q']<.05);basecorrect=(calls['v16']['q']<.05)&(np.sign(calls['v16']['lfc'])==np.sign(truth['lfc']));chosen='confidence' if 'confidence' in calls else names[1];newcorrect=(calls[chosen]['q']<.05)&(np.sign(calls[chosen]['lfc'])==np.sign(truth['lfc']))
    diagnostic={'family':family,'target':target,'seed':seed,'eligible_covered_genes':int(mask.sum()),'source_expert_error_correlation':cor,'baseline_raw_effect_L1_error':float(abs(error0).mean()),'PIE_raw_effect_L1_error':float(abs(error1).mean()),'confidence_changed_coordinates':int(np.count_nonzero(effects['confidence']!=effect)),'corrects_true_DE':int((true&~basecorrect&newcorrect).sum()),'breaks_true_DE':int((true&basecorrect&~newcorrect).sum()),'new_wrong_significant':int((mask&(calls[chosen]['q']<.05)&(np.sign(calls[chosen]['lfc'])!=np.sign(truth['lfc']))).sum())}
    prov={'family':family,'target':target,'seed':seed,'held_family':family,'actual_unique_cells':actual.shape[0],'template_unique_cells':400,'compared_populations':len(names)};np.savez_compressed(file,**arrays);write(side,{'inputs':inputs,'sha256':b.sha256(file),'records':local,'provenance':prov,'diagnostic':diagnostic});print('V18',stage,family,target,seed,'completed',flush=True);return clean(local),arrays,prov,diagnostic
   with ThreadPoolExecutor(max_workers=workers) as pool:outputs=list(pool.map(one,cells))
   for rr,_,prov,diag in outputs:scores.extend(rr);provenance.append(prov);diagnostics.append(diag)
   tp=np.asarray([a['truth'] for _,a,_,_ in outputs]);tj=np.asarray([a['truth_jack'] for _,a,_,_ in outputs]);axis=~np.isin(genes,targets);np.testing.assert_allclose(pds(np.broadcast_to(cp,tp.shape)[:,axis],tp[:,axis],cp[axis]),.5,atol=1e-12)
   for name in names:
    pp=np.asarray([a[name+'_pred'] for _,a,_,_ in outputs]);pj=np.asarray([a[name+'_jack'] for _,a,_,_ in outputs]);metric=panel_metrics(pp,tp,cp,pj,tj,cj,genes,[c['target'] for c in cells],targets);f=pd.DataFrame([r for rr,_,_,_ in outputs for r in rr if r['candidate']==name]);panels.append({**f.mean(numeric_only=True).to_dict(),'family':family,'context':ctx['context'],'seed':seed,'candidate':name,**metric})
   pd.DataFrame(scores).to_csv(A/f'{stage}_scores.csv',index=False);pd.DataFrame(panels).to_csv(A/f'{stage}_panels.csv',index=False);pd.DataFrame(diagnostics).to_csv(A/f'{stage}_expert_diagnostics.csv',index=False)
  del x,population;gc.collect()
 frame=pd.DataFrame(panels);summary=frame.groupby(['candidate','family'])[FIELDS].mean().groupby('candidate').mean();summary.to_csv(A/f'{stage}_summary.csv');write(A/f'{stage}_evaluation.json',{'complete':True,'population_comparisons':len(scores),'target_seed_cases':len(provenance),'provenance':provenance,'scores_sha256':b.sha256(A/f'{stage}_scores.csv'),'panels_sha256':b.sha256(A/f'{stage}_panels.csv'),'plan_sha256':b.sha256(A/'plan.json')});print(summary.to_string(),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);p.add_argument('--workers',type=int,default=4);a=p.parse_args();assert 1<=a.workers<=6
 with threadpool_limits(limits=2):evaluate(a.stage,a.workers)
