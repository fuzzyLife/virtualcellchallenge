"""V17 reserved verification with corrected, frozen disjoint target identities."""
import argparse,json,gc,hashlib
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
from scipy.stats import wasserstein_distance
from threadpoolctl import threadpool_limits
import build_submission as b
from v17_methods import ROOT,OUT,A,RAW_FAMILIES,CANDIDATES,fit_source,pool_variation,emit,strengths,logcpm
from fit_v7_response import read_cache,raw_inputs
from v16_methods import source_cache
from v8_methods import gated,panel_metrics
from v6_methods import moments
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked,pds
from evaluate_atlas_transfer import load_atlas_model
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio','shape_wasserstein']

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
 plan=json.loads((A/'confirmation_plan.json').read_text())
 assert plan['parent_plan_sha256']==b.sha256(A/'plan.json')
 for n,h in plan['dependencies'].items():assert b.sha256(ROOT/n)==h,'Frozen V17 dependency changed: '+n
 assert b.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')==plan['source_cache_sha256']
 return plan

def fit():
 plan=load_plan();data=read_cache();folder=A/'models';folder.mkdir(exist_ok=True)
 for f in RAW_FAMILIES:
  if not (A/f'{f}_variation_fits.csv').exists():fit_source(f,data['genes'],data['targets'])
 sources={f:b.sha256(A/f'{f}_variation_fits.csv') for f in RAW_FAMILIES}
 for f in RAW_FAMILIES+['all']:
  excluded=[] if f=='all' else [f];model=pool_variation(excluded,data['targets']);path=folder/(f+'.npz')
  if path.exists():raise FileExistsError('Preserve fitted V17 model')
  np.savez_compressed(path,**model,targets=data['targets']);write(path.with_suffix('.json'),{'excluded_families':excluded,'training_families':[v for v in RAW_FAMILIES if v not in excluded],'training_sources':{k:v for k,v in sources.items() if k not in excluded},'sha256':b.sha256(path),'plan_sha256':b.sha256(A/'plan.json'),'supported_targets':int((model['support']>0).sum()),'nonzero_targets':int((model['cv2']>0).sum()),'mean_cv2':float(model['cv2'].mean()),'max_cv2':float(model['cv2'].max())})
  print('V17 pooled',f,'nonzero',int((model['cv2']>0).sum()),'CV2 mean/max',float(model['cv2'].mean()),float(model['cv2'].max()),flush=True)
 write(A/'models_complete.json',{'complete':True,'plan_sha256':b.sha256(A/'plan.json'),'sources':sources,'models':{f:b.sha256(folder/(f+'.npz')) for f in RAW_FAMILIES+['all']}})

def load_model(family):
 file=A/'models'/(family+'.npz');rec=json.loads(file.with_suffix('.json').read_text());assert rec['plan_sha256']==b.sha256(A/'plan.json') and rec['sha256']==b.sha256(file)
 for f,h in rec['training_sources'].items():assert b.sha256(A/f'{f}_variation_fits.csv')==h
 assert family=='all' or family not in rec['training_families']
 with np.load(file,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def evaluate(stage,workers):
 plan=load_plan();data=read_cache();genes=data['genes'];targets=data['targets'];atlas,_=load_atlas_model();scores=[];panels=[];provenance=[]
 names=CANDIDATES if stage=='development' else (['v16'] if stage=='null' else ['v16',json.loads((A/'selection.json').read_text())['candidate'],'generic','shuffled'])
 cache=A/(stage+'_cache');cache.mkdir(exist_ok=True)
 for ctx in plan['contexts']:
  family=ctx['family'];x,meta,measured,_=raw_inputs(family,genes);population=x[ctx['left']];_,p,bulk,_=moments(population);query=p*1e6;query[~measured]=np.nan;ref=x[ctx['right']];cp,cj=profile_jackknife_chunked(ref)
  bank=source_cache(data,[family]);responses=np.asarray([gated(v.copy(),np.nan_to_num(query),int(o)) for v,o in zip(bank['base'],bank['own'])],np.float32);model=load_model(family)
  fitlog=logcpm(population);mu=np.asarray(fitlog.mean(0)).ravel();var=np.maximum(np.asarray(fitlog.power(2).mean(0)).ravel()-mu*mu,0);del fitlog
  byseed={}
  for c in (ctx['null_cases'] if stage=='null' else [r for r in ctx['cells'] if r['split']==stage]):byseed.setdefault(c['seed'],[]).append(c)
  for seed,cells in byseed.items():
   def one(cell):
    null=stage=='null';target='Non-Targeting' if null else cell['target'];ti=0 if null else int(np.flatnonzero(targets==target)[0]);own=None if null else int(bank['own'][ti]);effect=np.zeros(len(genes)) if null else responses[ti];file=cache/f'{family}_{target}_{seed}.npz';side=file.with_suffix('.json');inputs={'plan':b.sha256(A/'confirmation_plan.json'),'model':b.sha256(A/'models'/(family+'.npz')),'names':names}
    if side.exists():
     rec=json.loads(side.read_text());assert rec['inputs']==inputs and rec['sha256']==b.sha256(file)
     with np.load(file,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
     return rec['records'],arrays,rec['provenance']
    actual=x[cell['actual']];template=x[cell['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64);assert len(set(cell['template']))==400 and not set(cell['actual'])&set(cell['template']);assert not set(cell['template'])&set(ctx['right'])
    truth=de(actual,ref,genes);tp,tj=profile_jackknife_chunked(actual);arrays={'truth':tp,'truth_jack':tj};local=[]
    selection=np.flatnonzero(measured&~np.isin(genes,targets)&(mu>.25)&(var>.01));selection=selection[np.argsort(abs(effect[selection])/np.sqrt(var[selection]+.1))[-256:]];w=effect[selection]/(var[selection]+.1);w/=max(np.linalg.norm(w),1e-12);actual_projection=np.asarray(logcpm(actual)[:,selection]@w).ravel()
    for name in names:
     cv=strengths(model,name,ti);count,diag=emit(template,np.arange(400),effect,own,seed+ti,cv,prepared=(libs,p,bulk),atlas=atlas)
     if null:
      other,_=emit(template,np.arange(400),effect,None,seed+ti,.25,prepared=(libs,p,bulk),atlas=atlas);assert (count!=other).nnz==0
     pp,pj=profile_jackknife_chunked(count);arrays[name+'_pred']=pp;arrays[name+'_jack']=pj;metrics=direction_metrics(de(count,ref,genes),truth,genes,-1 if own is None else own)
     pred_projection=np.asarray(logcpm(count)[:,selection]@w).ravel();shape=float(wasserstein_distance(actual_projection,pred_projection))
     local.append({'family':family,'context':ctx['context'],'target':target,'seed':seed,'candidate':name,'actual_cells':actual.shape[0],'prediction_cells':400,'reference_cells':ref.shape[0],'count_sha256':digest(count),'shape_wasserstein':shape,'actual_projection_variance':float(actual_projection.var()),'prediction_projection_variance':float(pred_projection.var()),'learned_cv2':cv,'null_bitwise_v16_identity':bool(null),**metrics,**diag})
    prov={'family':family,'target':target,'seed':seed,'held_family':family,'actual_unique_cells':actual.shape[0],'template_unique_cells':400,'compared_populations':len(names)};np.savez_compressed(file,**arrays);write(side,{'inputs':inputs,'sha256':b.sha256(file),'records':local,'provenance':prov});print('V17',stage,family,target,seed,'completed',flush=True);return clean(local),arrays,prov
   with ThreadPoolExecutor(max_workers=workers) as pool:outputs=list(pool.map(one,cells))
   for rr,_,prov in outputs:scores.extend(rr);provenance.append(prov)
   tp=np.asarray([a['truth'] for _,a,_ in outputs]);tj=np.asarray([a['truth_jack'] for _,a,_ in outputs])
   if stage!='null':
    axis=~np.isin(genes,targets);np.testing.assert_allclose(pds(np.broadcast_to(cp,tp.shape)[:,axis],tp[:,axis],cp[axis]),.5,atol=1e-12)
   for name in names:
    pp=np.asarray([a[name+'_pred'] for _,a,_ in outputs]);pj=np.asarray([a[name+'_jack'] for _,a,_ in outputs]);metric={} if stage=='null' else panel_metrics(pp,tp,cp,pj,tj,cj,genes,[c['target'] for c in cells],targets);f=pd.DataFrame([r for rr,_,_ in outputs for r in rr if r['candidate']==name]);panels.append({**f.mean(numeric_only=True).to_dict(),'family':family,'context':ctx['context'],'seed':seed,'candidate':name,**metric})
   pd.DataFrame(scores).to_csv(A/f'{stage}_scores.csv',index=False);pd.DataFrame(panels).to_csv(A/f'{stage}_panels.csv',index=False)
  del x,population;gc.collect()
 frame=pd.DataFrame(panels);fields=['predicted_calls','reference_calls','shape_wasserstein'] if stage=='null' else FIELDS
 summary=frame.groupby(['candidate','family'])[fields].mean().groupby('candidate').mean();summary.to_csv(A/f'{stage}_summary.csv');write(A/f'{stage}_evaluation.json',{'complete':True,'population_comparisons':len(scores),'target_seed_cases':len(provenance),'provenance':provenance,'scores_sha256':b.sha256(A/f'{stage}_scores.csv'),'panels_sha256':b.sha256(A/f'{stage}_panels.csv'),'plan_sha256':b.sha256(A/'confirmation_plan.json')});print(summary.to_string(),flush=True)
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['confirmation']);parser.add_argument('--workers',type=int,default=6);args=parser.parse_args()
 if not 1<=args.workers<=6:parser.error('workers must be 1..6')
 with threadpool_limits(limits=2):
  if args.stage=='models':fit()
  else:evaluate(args.stage,args.workers)
