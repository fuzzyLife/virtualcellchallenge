"""Execute V22 oracle, conditional fitting and paired generated-count tests."""
import argparse,json,gc,time,sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy import sparse
from scipy.stats import wasserstein_distance
from threadpoolctl import threadpool_limits
from v22_methods import ROOT,A,FAMILIES,NAMES,load_plan,fit_model,load_model,raw,emit,write,digest,logcpm,query_states
from fit_v7_response import read_cache
from v16_methods import source_arrays
from v6_methods import moments
from v8_methods import gated,panel_metrics
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from evaluate_atlas_transfer import load_atlas_model
import build_submission as b
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio','shape_wasserstein']
ORACLE=['pred_mean_compact','measured_mean_compact','pred_mean_treated_shape','measured_mean_treated_shape','real_replica']

def evaluate(stage,workers):
 plan=load_plan();data=read_cache();genes=data['genes'];targets=data['targets'];atlas,_=load_atlas_model();records=[];panels=[];contexts=plan['contexts']
 if stage=='confirmation':
  selected=json.loads((A/'selection.json').read_text());assert selected['new_candidate_passes']
  contexts=[json.loads((A/'reserve_plan.json').read_text())]
 names=ORACLE if stage=='oracle' else (['v16','conditional'] if stage=='null' else NAMES)
 for ctx in contexts:
  family=ctx['family'];x,meta,measured,source=raw(family,genes);population=x[ctx['left']];_,p,bulk,_=moments(population);ref=x[ctx['right']];cp,cj=profile_jackknife_chunked(ref);bank=source_arrays(data,[family]);effects=np.array([gated(v.copy(),p*1e6,int(o)) for v,o in zip(bank['base'],bank['own'])],np.float32)
  model=None if stage=='oracle' else load_model(family);zs=None if stage=='oracle' else query_states(population,model);lookup={int(v):i for i,v in enumerate(ctx['left'])}
  if stage=='oracle':cells=ctx['oracle']
  elif stage=='null':
   ids=np.random.default_rng(22090+FAMILIES.index(family)).permutation(ctx['left']);cells=[{'target':'Non-Targeting','seed':22090,'template':ids[:400].tolist(),'actual':ids[400:600].tolist()}]
  else:cells=ctx['cells']
  group={}
  for c in cells:group.setdefault(c['seed'] if stage!='oracle' else 0,[]).append(c)
  cache=A/(stage+'_cache');cache.mkdir(exist_ok=True)
  for seed,subset in group.items():
   def one(cell):
    target=cell['target'];null=stage=='null';ti=0 if null else int(np.flatnonzero(targets==target)[0]);own=None if null else int(bank['own'][ti]);effect=np.zeros(len(genes)) if null else effects[ti]
    file=cache/f'{family}_{target}_{cell["seed"]}.npz';side=file.with_suffix('.json');inputs={'plan':b.sha256(A/'plan.json'),'names':names,'model':None if model is None else b.sha256(A/'models'/(family+'.npz')),'reserve':b.sha256(A/'reserve_plan.json') if stage=='confirmation' else None}
    if side.exists():
     rec=json.loads(side.read_text());assert rec['inputs']==inputs and rec['sha256']==b.sha256(file)
     with np.load(file,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
     return rec['records'],arrays
    actual=x[cell['actual']];template=x[cell['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64)
    assert len(set(cell['actual']))==len(cell['actual']) and len(set(cell['template']))==len(cell['template']);assert not set(cell['actual'])&set(cell['template']);assert not (set(cell['actual'])|set(cell['template']))&set(ctx['right'])
    truth=de(actual,ref,genes);tp,tj=profile_jackknife_chunked(actual);arrays={'truth':tp,'truth_jack':tj};local=[]
    axis=np.flatnonzero(measured&~np.isin(genes,targets)&(p*1e6>5));axis=axis[np.argsort(abs(effect[axis]),kind='stable')[-256:]];direction=effect[axis]/max(np.linalg.norm(effect[axis]),1e-9);trueproj=np.asarray(logcpm(actual)[:,axis]@direction).ravel();trueproj-=trueproj.mean()
    calibration=None
    if stage=='oracle':
     assert not set(cell['calibration'])&set(cell['actual']);calibration=x[cell['calibration']];libs,op,ob,_=moments(calibration);libs=libs.astype(np.int64);factor=np.exp2(effect);factor[own]=.1;pp=p*factor;pp/=pp.sum();bb=bulk*factor;bb/=bb.sum();cprob=calibration.toarray()/libs[:,None]
    z=None if model is None else zs[[lookup[int(v)] for v in cell['template']]]
    for name in names:
     if stage=='oracle':
      if name=='real_replica':count=calibration.copy();diag={}
      else:
       wp,wb=(op,ob) if name.startswith('measured') else (pp,bb)
       if name.endswith('treated_shape'):initial=cprob*np.divide(wp,op,out=np.ones_like(wp),where=op>0)[None,:]
       else:
        expected=libs[:,None]*wp;empirical=template.toarray().astype(float)*np.divide(wp,p,out=np.zeros_like(wp),where=p>0)[None,:];empirical*=np.divide(libs,np.maximum(empirical.sum(1),1))[:,None];initial=.25*expected+.65*np.random.default_rng(cell['seed']).poisson(expected)+.1*empirical
       initial/=np.maximum(initial.sum(1,keepdims=True),1e-12);count=sparse.csr_matrix(atlas.dual_moment_counts(initial,wp,wb,depths=libs,seed=cell['seed'],iterations=200));diag={'mean_probability_L1':float(abs(moments(count)[1]-wp).sum()),'bulk_probability_L1':float(abs(moments(count)[2]-wb).sum())}
     else:count,diag=emit(template,effect,own,cell['seed']+ti,name,(libs,p,bulk),atlas,model,ti,z)
     np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),libs);assert np.issubdtype(count.dtype,np.integer) and (count.data>=0).all()
     pp1,pj=profile_jackknife_chunked(count);arrays[name+'_pred']=pp1;arrays[name+'_jack']=pj;metrics=direction_metrics(de(count,ref,genes),truth,genes,-1 if own is None else own);proj=np.asarray(logcpm(count)[:,axis]@direction).ravel();proj-=proj.mean();shape=float(wasserstein_distance(trueproj,proj))
     local.append({'family':family,'target':target,'seed':cell['seed'],'candidate':name,'actual_cells':actual.shape[0],'prediction_cells':count.shape[0],'reference_cells':ref.shape[0],'count_sha256':digest(count),'shape_wasserstein':shape,**metrics,**diag})
     if null and name=='v16':base=count.copy()
     if null and name=='conditional':assert (base!=count).nnz==0
    np.savez_compressed(file,**arrays);write(side,{'inputs':inputs,'sha256':b.sha256(file),'records':local,'original_input_record':source});print('V22',stage,family,target,cell['seed'],'complete',flush=True);return local,arrays
   with ThreadPoolExecutor(max_workers=workers) as pool:outputs=list(pool.map(one,subset))
   for rr,_ in outputs:records.extend(rr)
   tp=np.array([a['truth'] for _,a in outputs]);tj=np.array([a['truth_jack'] for _,a in outputs])
   for name in names:
    pp=np.array([a[name+'_pred'] for _,a in outputs]);pj=np.array([a[name+'_jack'] for _,a in outputs]);m={} if stage=='null' else panel_metrics(pp,tp,cp,pj,tj,cj,genes,[c['target'] for c in subset],targets);frame=pd.DataFrame([r for rr,_ in outputs for r in rr if r['candidate']==name]);panels.append({'family':family,'seed':seed,'candidate':name,**frame.mean(numeric_only=True).to_dict(),**m})
   pd.DataFrame(records).to_csv(A/f'{stage}_scores.csv',index=False);pd.DataFrame(panels).to_csv(A/f'{stage}_panels.csv',index=False)
  del x,population;gc.collect()
 frame=pd.DataFrame(panels);fields=['predicted_calls','reference_calls'] if stage=='null' else FIELDS;summary=frame.groupby(['candidate','family'])[fields].mean().groupby('candidate').mean();summary.to_csv(A/f'{stage}_summary.csv');write(A/f'{stage}_evaluation.json',{'complete':True,'populations':len(records),'case_count':len(records)//len(names),'plan_sha256':b.sha256(A/'plan.json'),'scores_sha256':b.sha256(A/f'{stage}_scores.csv'),'panels_sha256':b.sha256(A/f'{stage}_panels.csv'),'scope':'Held-family generated-count test; public raw metrics, missing values retained; oracle diagnostic uses held treatments and is not deployable.'});print(summary.to_string(),flush=True)

def main():
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['oracle','models','null','development','confirmation']);p.add_argument('--workers',type=int,default=4);args=p.parse_args();assert 1<=args.workers<=6
 with threadpool_limits(limits=2):
  if args.stage=='models':
   for f in FAMILIES+['CD4','all']:fit_model(f)
  else:evaluate(args.stage,args.workers)
if __name__=='__main__':main()
