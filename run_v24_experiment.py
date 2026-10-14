"""Paired V24 count tests; fresh condition only after locking development choice."""
from pathlib import Path
import argparse,json,gc,hashlib,shutil
from concurrent.futures import ProcessPoolExecutor,ThreadPoolExecutor,as_completed
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from threadpoolctl import threadpool_limits
from scipy.stats import wasserstein_distance
from v24_methods import ROOT,A,FAMILIES,NAMES,load_plan,make_teacher,raw,predict,emit,write,digest,logcpm
from fit_v7_response import read_cache
from v6_methods import moments
from v8_methods import panel_metrics
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from evaluate_atlas_transfer import load_atlas_model
import build_submission as b
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio','shape_wasserstein']

def panel(task):
 ctx,seed,stage=task;family=ctx['family'];plan=load_plan();data=read_cache();genes=data['genes'];targets=data['targets'];teacher=make_teacher(False);x,meta,measured,source=raw(family,genes);population=x[ctx['left']];_,p,bulk,_=moments(population);query=p*1e6;query[~measured]=np.nan;effects,diagnostic=predict(data,teacher,query,ctx['excluded']);ref=x[ctx['right']];cp,cj=profile_jackknife_chunked(ref);atlas,_=load_atlas_model();cells=[v for v in ctx['cells'] if v['seed']==seed];folder=A/'parallel'/f'{family}_{seed}_{stage}';folder.mkdir(parents=True,exist_ok=True);cache=folder/'count_cache';cache.mkdir(exist_ok=True)
 def one(cell):
  target=cell['target'];ti=int(np.flatnonzero(targets==target)[0]);own=int(teacher['own'][ti]);file=cache/f'{family}_{target}_{seed}.npz';side=file.with_suffix('.json');inputs={'plan':b.sha256(A/'plan.json'),'teacher':b.sha256(A/'models/teacher_development.npz'),'reserve':b.sha256(A/'reserve_plan.json') if stage=='confirmation' else None,'names':NAMES}
  if side.exists():
   rec=json.loads(side.read_text());assert rec['inputs']==inputs and rec['sha256']==b.sha256(file)
   with np.load(file,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
   return rec['records'],arrays
  actual=x[cell['actual']];template=x[cell['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64);assert len(set(cell['actual']))==actual.shape[0] and len(set(cell['template']))==template.shape[0];assert not set(cell['actual'])&(set(ctx['left'])|set(ctx['right']));assert set(cell['template'])<=set(ctx['left']);truth=de(actual,ref,genes);tp,tj=profile_jackknife_chunked(actual);arrays={'truth':tp,'truth_jack':tj};local=[];by_effect={};axis=np.flatnonzero(measured&~np.isin(genes,targets)&(p*1e6>5));axis=axis[np.argsort(abs(effects['v16'][ti,axis]),kind='stable')[-256:]];direction=effects['v16'][ti,axis]/max(np.linalg.norm(effects['v16'][ti,axis]),1e-9);projection=np.asarray(logcpm(actual)[:,axis]@direction).ravel();projection-=projection.mean()
  for name in NAMES:
   effect=effects[name][ti];key=np.where(effect==0,0.,effect).astype(np.float32).tobytes()
   if key in by_effect:
    row,pp,pj=by_effect[key];row=dict(row,candidate=name,reused_equal_effect=True)
   else:
    count,diag=emit(template,effect,own,seed+ti,(libs,p,bulk),atlas);np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),libs);pp,pj=profile_jackknife_chunked(count);metrics=direction_metrics(de(count,ref,genes),truth,genes,own);z=np.asarray(logcpm(count)[:,axis]@direction).ravel();z-=z.mean();row={'family':family,'target':target,'seed':seed,'candidate':name,'actual_cells':actual.shape[0],'prediction_cells':count.shape[0],'reference_cells':ref.shape[0],'count_sha256':digest(count),'reused_equal_effect':False,'shape_wasserstein':float(wasserstein_distance(projection,z)),**metrics,**diag};by_effect[key]=(row,pp,pj)
   arrays[name+'_pred']=pp;arrays[name+'_jack']=pj;local.append(row)
  np.savez_compressed(file,**arrays);write(side,{'inputs':inputs,'sha256':b.sha256(file),'records':local,'source':source,'query_routing':diagnostic});print('V24',stage,family,target,seed,'complete',flush=True);return local,arrays
 with threadpool_limits(limits=2),ThreadPoolExecutor(max_workers=4) as pool:outputs=list(pool.map(one,cells))
 scores=[r for rows,_ in outputs for r in rows];tp=np.array([z['truth'] for _,z in outputs]);tj=np.array([z['truth_jack'] for _,z in outputs]);panels=[]
 for name in NAMES:
  pp=np.array([z[name+'_pred'] for _,z in outputs]);pj=np.array([z[name+'_jack'] for _,z in outputs]);metrics=panel_metrics(pp,tp,cp,pj,tj,cj,genes,[c['target'] for c in cells],targets);f=pd.DataFrame([r for r in scores if r['candidate']==name]);panels.append({'family':family,'seed':seed,'candidate':name,**f.mean(numeric_only=True).to_dict(),**metrics})
 pd.DataFrame(scores).to_csv(folder/'scores.csv',index=False);pd.DataFrame(panels).to_csv(folder/'panels.csv',index=False);write(folder/'complete.json',{'complete':True,'query_routing':diagnostic,'cases':len(cells),'scores_sha256':b.sha256(folder/'scores.csv'),'panels_sha256':b.sha256(folder/'panels.csv')});return str(folder)

def run(stage):
 plan=load_plan();make_teacher(False)
 if stage=='confirmation':
  selected=json.loads((A/'selection.json').read_text());assert selected['new_candidate_passes'];start=A/'reserve_execution_started.json'
  if start.exists():raise FileExistsError('Reserve already opened; inspect identical caches before any operational recovery')
  write(start,{'selection_sha256':b.sha256(A/'selection.json'),'reserve_plan_sha256':b.sha256(A/'reserve_plan.json'),'runner_sha256':b.sha256(Path(__file__))});contexts=[json.loads((A/'reserve_plan.json').read_text())]
 else:contexts=plan['contexts']
 tasks=[(ctx,seed,stage) for ctx in contexts for seed in sorted({c['seed'] for c in ctx['cells']})];folders=[]
 with ProcessPoolExecutor(max_workers=min(6,len(tasks))) as pool:
  for future in as_completed([pool.submit(panel,t) for t in tasks]):folders.append(Path(future.result()))
 scores=pd.concat([pd.read_csv(f/'scores.csv') for f in folders],ignore_index=True).sort_values(['family','seed','target','candidate']);panels=pd.concat([pd.read_csv(f/'panels.csv') for f in folders],ignore_index=True).sort_values(['family','seed','candidate']);scores.to_csv(A/f'{stage}_scores.csv',index=False);panels.to_csv(A/f'{stage}_panels.csv',index=False);summary=panels.groupby(['candidate','family'])[FIELDS].mean().groupby('candidate').mean();summary.to_csv(A/f'{stage}_summary.csv');write(A/f'{stage}_evaluation.json',{'complete':True,'populations':len(scores),'cases':len(scores)//len(NAMES),'plan_sha256':b.sha256(A/'plan.json'),'scores_sha256':b.sha256(A/f'{stage}_scores.csv'),'panels_sha256':b.sha256(A/f'{stage}_panels.csv'),'source_contexts':[c['family'] for c in contexts],'runner_sha256':b.sha256(Path(__file__))});print(summary.to_string(),flush=True)

def main():
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['models','development','confirmation']);args=p.parse_args()
 with threadpool_limits(limits=2):
  if args.stage=='models':make_teacher(False)
  else:run(args.stage)
if __name__=='__main__':main()
