"""Paired V26 count calibration and treatment tests without changing response means."""
from pathlib import Path
import argparse,json,gc
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy.stats import wasserstein_distance
from threadpoolctl import threadpool_limits
from v26_methods import ROOT,A,NAMES,KINDS,family_key,load_plan,load_decoder,raw,baseline,moments,state_values,descriptor,emit,digest,logcpm,write
from fit_v7_response import read_cache
from v12_cpu_metrics import de,direction_metrics
from v8_methods import panel_metrics
from audit_transfer_validation import profile_jackknife_chunked
from evaluate_atlas_transfer import load_atlas_model
import build_submission as b
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio','shape_wasserstein']

def distribution(x,ref,axis):
 y=logcpm(x)[:,axis];r=logcpm(ref)[:,axis];dx=np.asarray((y>0).mean(0)).ravel();dr=np.asarray((r>0).mean(0)).ravel();vx=np.maximum(np.asarray(y.power(2).mean(0)).ravel()-np.asarray(y.mean(0)).ravel()**2,0);vr=np.maximum(np.asarray(r.power(2).mean(0)).ravel()-np.asarray(r.mean(0)).ravel()**2,0)
 return {'detection_MAE':float(abs(dx-dr).mean()),'log_variance_MAE':float(abs(np.log(vx+.01)-np.log(vr+.01)).mean())}

def main(stage,family=None):
 from v24_runtime import require_science
 require_science();plan=load_plan();data=read_cache();genes=data['genes'];targets=data['targets'];atlas,_=load_atlas_model();contexts=[c for c in plan['contexts'] if family is None or family_key(c['family'])==family];records=[];panels=[];tag='null_'+family if stage=='null' else stage;assert not (A/f'{tag}_evaluation.json').exists()
 for ctx in contexts:
  key=family_key(ctx['family']);models={k:load_decoder(key,k) for k in KINDS};x,meta,measured,source=raw(ctx['family'],genes);ids=ctx['null']['fit_controls'] if stage=='null' else ctx['left'];population=x[ids];_,p,bulk,_=moments(population);query=p*1e6;query[~measured]=np.nan;effects=baseline(data,query,ctx['excluded']);own=np.array([np.flatnonzero(genes==t)[0] for t in targets]);refids=ctx['null']['reference'] if stage=='null' else ctx['right'];ref=x[refids];cp,cj=profile_jackknife_chunked(ref);states=state_values(x,meta,ids,models['conditional']['state_axis'],models['conditional']['state_basis'],models['conditional']['state_scale']);lookup={int(v):i for i,v in enumerate(ids)};perm=np.random.default_rng(26026).permutation(300);cells=[ctx['null']] if stage=='null' else ctx['cells'];cache=A/'parallel'/tag/ctx['family'];cache.mkdir(parents=True,exist_ok=True)
  def one(cell):
   target=cell['target'];null=stage=='null';ti=0 if null else int(np.flatnonzero(targets==target)[0]);o=None if null else int(own[ti]);effect=np.zeros(len(genes)) if null else effects[ti];template=x[cell['template']];actual=x[cell['actual']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64);state=states[[lookup[int(v)] for v in cell['template']]];side=cache/f'{target}_{cell["seed"]}.json';file=side.with_suffix('.npz');inputs={'plan':b.sha256(A/'plan.json'),'models':{k:models[k]['record']['sha256'] for k in KINDS}}
   if side.exists():
    rec=json.loads(side.read_text());assert rec['inputs']==inputs and rec['sha256']==b.sha256(file)
    with np.load(file,allow_pickle=False) as z:return rec['records'],{n:z[n].copy() for n in z.files}
   assert len(set(cell['actual']))==len(cell['actual']) and len(set(cell['template']))==len(cell['template']);assert not set(cell['actual'])&set(ids) and not set(cell['actual'])&set(refids);assert set(cell['template'])<=set(ids)
   truth=de(actual,ref,genes);tp,tj=profile_jackknife_chunked(actual);arrays={'truth':tp,'truth_jack':tj};axis=np.flatnonzero(measured&~np.isin(genes,targets)&(p*1e6>5));project=axis[np.argsort(abs(effect[axis]),kind='stable')[-256:]];direction=effect[project]/max(np.linalg.norm(effect[project]),1e-9);tproj=np.asarray(logcpm(actual)[:,project]@direction).ravel();tproj-=tproj.mean();local=[];reference=None;conditional_count=None;conditional_diag=None
   for name in NAMES:
    model=models['conditional' if name in ['v24','empirical','target_shuffled'] else name];condition=None if name!='target_shuffled' or null else descriptor(effects[perm[ti]],p,model['geometry']);seed=cell['seed']+ti
    if null and name=='target_shuffled':count=conditional_count.copy();diag=dict(conditional_diag)
    else:count,diag=emit(template,effect,o,seed,name,(libs,p,bulk),atlas,model,state,condition,reference)
    if name=='v24':reference=count
    if name=='conditional':conditional_count=count;conditional_diag=diag
    np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),libs);assert np.issubdtype(count.dtype,np.integer) and (count.data>=0).all();pp,pj=profile_jackknife_chunked(count);arrays[name+'_pred']=pp;arrays[name+'_jack']=pj;d=de(count,ref,genes);metrics=direction_metrics(d,truth,genes,-1 if o is None else o);proj=np.asarray(logcpm(count)[:,project]@direction).ravel();proj-=proj.mean();local.append({'family':ctx['family'],'target':target,'seed':cell['seed'],'candidate':name,'actual_cells':actual.shape[0],'prediction_cells':count.shape[0],'reference_cells':ref.shape[0],'count_sha256':digest(count),'shape_wasserstein':float(wasserstein_distance(tproj,proj)),'null_DE_calls':int(((d['q']<.05)&d['eligible']).sum()),'real_null_DE_calls':int(((truth['q']<.05)&truth['eligible']).sum()),**distribution(count,ref,axis),**{'real_'+k:v for k,v in distribution(actual,ref,axis).items()},**metrics,**diag})
    print('V26',tag,ctx['family'],target,name,'count complete',flush=True)
   np.savez_compressed(file,**arrays);write(side,{'inputs':inputs,'sha256':b.sha256(file),'records':local,'source':source});return local,arrays
  with ThreadPoolExecutor(max_workers=2 if stage!='null' else 1) as pool:outputs=list(pool.map(one,cells))
  records.extend(r for rows,_ in outputs for r in rows)
  for seed in sorted({c['seed'] for c in cells}):
   use=[i for i,c in enumerate(cells) if c['seed']==seed];tp=np.array([outputs[i][1]['truth'] for i in use]);tj=np.array([outputs[i][1]['truth_jack'] for i in use])
   for name in NAMES:
    pp=np.array([outputs[i][1][name+'_pred'] for i in use]);pj=np.array([outputs[i][1][name+'_jack'] for i in use]);metrics={} if stage=='null' else panel_metrics(pp,tp,cp,pj,tj,cj,genes,[cells[i]['target'] for i in use],targets);f=pd.DataFrame([r for i in use for r in outputs[i][0] if r['candidate']==name]);panels.append({'family':ctx['family'],'seed':seed,'candidate':name,**f.mean(numeric_only=True).to_dict(),**metrics})
  del x,population;gc.collect()
 pd.DataFrame(records).to_csv(A/f'{tag}_scores.csv',index=False);pd.DataFrame(panels).to_csv(A/f'{tag}_panels.csv',index=False);write(A/f'{tag}_evaluation.json',{'complete':True,'cases':len(records)//len(NAMES),'populations':len(records),'plan_sha256':b.sha256(A/'plan.json'),'scores_sha256':b.sha256(A/f'{tag}_scores.csv'),'panels_sha256':b.sha256(A/f'{tag}_panels.csv')});print('V26',tag,'complete',len(records),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['null','development']);p.add_argument('--family');args=p.parse_args()
 with threadpool_limits(limits=2):main(args.stage,args.family)
