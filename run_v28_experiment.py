"""V28 sampling-stage ablations and paired control/treatment populations."""
import argparse,json,gc,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy import sparse
from scipy.stats import wasserstein_distance
from threadpoolctl import threadpool_limits
from v28_methods import ROOT,A,ARMS,arms_for,load_plan,load_model,emit,write,digest,logcpm,project
from v26_methods import raw,baseline,moments,state_values,descriptor,family_key,emit as absolute_emit,fit_probabilities,preserve_depth_margins,repair_cpm
from fit_v7_response import read_cache
from v12_cpu_metrics import de,direction_metrics
from run_v26_experiment import distribution
from v8_methods import panel_metrics
from audit_transfer_validation import profile_jackknife_chunked
from evaluate_atlas_transfer import load_atlas_model
import build_submission as b

def score_distribution(count,reference,actual,genes,axis):
 d=de(count,reference,genes);t=de(actual,reference,genes)
 return {'null_DE_calls':int(((d['q']<.05)&d['eligible']).sum()),'real_null_DE_calls':int(((t['q']<.05)&t['eligible']).sum()),**distribution(count,reference,axis),**{'real_'+k:v for k,v in distribution(actual,reference,axis).items()}}

def evaluate(stage,key):
 ARMS=arms_for(key)
 plan=load_plan();data=read_cache();genes=data['genes'];targets=data['targets'];atlas,_=load_atlas_model();models={'masked':load_model(key),'panel':load_model(key,True)};records=[];panels=[];tag=stage+'_'+key;contexts=[c for c in plan['contexts'] if family_key(c['family'])==key];perm=np.random.default_rng(27027).permutation(300)
 for ctx in contexts:
  x,meta,measured,source=raw(ctx['family'],genes);cells=ctx['null_cases'] if stage=='null' else ctx['cells'];folder=A/'parallel'/tag/ctx['family'];folder.mkdir(parents=True,exist_ok=True)
  def one(cell):
   null=stage=='null';target=cell['target'];ti=0 if null else int(np.flatnonzero(targets==target)[0]);seed=cell['seed']+ti;ids=cell['fit_controls'] if null else ctx['left'];refids=cell['reference'] if null else ctx['right'];fitids=cell['calibration_fit'] if null else [i for i in ids if i not in set(cell['template'])];side=folder/f'{target}_{cell["seed"]}.json';arrfile=side.with_suffix('.npz');inputhash={'plan':b.sha256(A/'plan.json'),'masked':models['masked']['record']['sha256'],'panel':models['panel']['record']['sha256']}
   if side.exists():
    rec=json.loads(side.read_text());assert rec['inputs']==inputhash and b.sha256(arrfile)==rec['sha256']
    with np.load(arrfile,allow_pickle=False) as z:return rec['records'],{n:z[n].copy() for n in z.files}
   assert len(set(cell['template']))==400;assert len(fitids)>=200
   ss=[set(cell['actual']),set(refids),set(cell['template']),set(fitids)]
   assert all(not ss[i]&ss[j] for i in range(4) for j in range(i));assert (meta.gene_target.to_numpy(str)[list(ss[1]|ss[2]|ss[3])]=='Non-Targeting').all()
   p,bulk=moments(x[ids])[1:3];query=p*1e6;query[~measured]=np.nan;effects=None if null else baseline(data,query,ctx['excluded']);effect=np.zeros(len(genes)) if null else effects[ti];own=None if null else int(np.flatnonzero(genes==target)[0]);template=x[cell['template']];actual=x[cell['actual']];ref=x[refids];fit=x[fitids];libs=template.sum(1).A1.astype(np.int64)
   state=state_values(x,meta,ids,models['masked']['state_axis'],models['masked']['state_basis'],models['masked']['state_scale']) if key=='H1' else np.zeros((len(ids),8),np.float32);lookup={int(v):i for i,v in enumerate(ids)};state=state[[lookup[int(v)] for v in cell['template']]]
   truth=de(actual,ref,genes);tp,tj=profile_jackknife_chunked(actual);arrays={'truth':tp,'truth_jack':tj};axis=np.flatnonzero(measured&~np.isin(genes,targets)&(p*1e6>5));pr=axis[np.argsort(abs(effect[axis]),kind='stable')[-256:]];direction=effect[pr]/max(np.linalg.norm(effect[pr]),1e-9);tproj=np.asarray(logcpm(actual)[:,pr]@direction).ravel();tproj-=tproj.mean();rows=[];base=None;counts={};diags={}
   for arm in ARMS:
    started=time.time();model=models['panel'] if arm in ['panel_relative','target_shuffled'] else models['masked'];override=descriptor(effects[perm[ti]],p,model['geometry']) if arm=='target_shuffled' and not null else None
    duplicate={'relative':'anchor_only','calibrated_relative':'calibrated_anchor','panel_relative':'calibrated_anchor','target_shuffled':'calibrated_anchor'}.get(arm) if null else None
    if duplicate:count=counts[duplicate];diag=diags[duplicate].copy()
    elif arm in ['v24','v26']:
     count,diag=absolute_emit(template,effect,own,seed,'v24' if arm=='v24' else 'conditional',(libs,p,bulk),atlas,model,state,reference=base)
    else:count,diag=emit(template,effect,own,p,seed,base,arm,model,state,fit,override)
    if arm=='v24':base=count
    counts[arm]=count;diags[arm]=diag
    assert count.shape==(400,len(genes)) and np.issubdtype(count.dtype,np.integer) and (count.data>=0).all();np.testing.assert_array_equal(count.sum(1).A1,libs)
    if arm!='v24':np.testing.assert_array_equal(count.sum(0).A1,base.sum(0).A1)
    pp,pj=profile_jackknife_chunked(count);arrays[arm+'_pred']=pp;arrays[arm+'_jack']=pj;d=de(count,ref,genes);metric=direction_metrics(d,truth,genes,-1 if own is None else own);proj=np.asarray(logcpm(count)[:,pr]@direction).ravel();proj-=proj.mean();row={'family':ctx['family'],'target':target,'seed':cell['seed'],'candidate':arm,'count_sha256':digest(count),'actual_cells':len(cell['actual']),'prediction_cells':400,'reference_cells':len(refids),'calibration_fit_cells':len(fitids),'shape_wasserstein':float(wasserstein_distance(tproj,proj)),**score_distribution(count,ref,actual,genes,axis),**metric,**diag};rows.append(row)
    if arm in ['v24','anchor_only','panel_relative','target_shuffled']:sparse.save_npz(folder/f'{target}_{cell["seed"]}_{arm}_counts.npz',count)
    print('V28',tag,ctx['family'],target,cell['seed'],arm,'calls',row['null_DE_calls'],'fidelity',row['fidelity'],'seconds',round(time.time()-started,1),flush=True)
   # Record actual metadata overlap; cell-disjoint is not guide/batch-disjoint.
   metadata={}
   for col in ['sample','guide_id','guide_target']:
    if col in meta:
     sets=[set(meta.iloc[v][col].astype(str)) for v in [fitids,cell['template'],cell['actual'],refids]];metadata[col]={'fit_unique':len(sets[0]),'template_unique':len(sets[1]),'actual_unique':len(sets[2]),'reference_unique':len(sets[3]),'fit_actual_overlap':len(sets[0]&sets[2]),'template_actual_overlap':len(sets[1]&sets[2])}
   np.savez_compressed(arrfile,**arrays);write(side,{'inputs':inputhash,'sha256':b.sha256(arrfile),'records':rows,'metadata_overlap':metadata,'calibration_fit_rows':fitids,'source':source});return rows,arrays
  with ThreadPoolExecutor(max_workers=2) as pool:outputs=list(pool.map(one,cells))
  records.extend(r for rows,_ in outputs for r in rows)
  for seed in sorted({c['seed'] for c in cells}):
   use=[i for i,c in enumerate(cells) if c['seed']==seed];tp=np.array([outputs[i][1]['truth'] for i in use]);tj=np.array([outputs[i][1]['truth_jack'] for i in use]);refids=cells[use[0]]['reference'] if stage=='null' else ctx['right'];cp,cj=profile_jackknife_chunked(x[refids])
   for arm in ARMS:
    pp=np.array([outputs[i][1][arm+'_pred'] for i in use]);pj=np.array([outputs[i][1][arm+'_jack'] for i in use]);metrics={} if stage=='null' else panel_metrics(pp,tp,cp,pj,tj,cj,genes,[cells[i]['target'] for i in use],targets);f=pd.DataFrame([r for i in use for r in outputs[i][0] if r['candidate']==arm]);panels.append({**f.mean(numeric_only=True).to_dict(),'family':ctx['family'],'seed':seed,'candidate':arm,**metrics})
  del x;gc.collect()
 pd.DataFrame(records).to_csv(A/f'{tag}_scores.csv',index=False);pd.DataFrame(panels).to_csv(A/f'{tag}_panels.csv',index=False);write(A/f'{tag}_evaluation.json',{'complete':True,'cases':len(records)//len(ARMS),'populations':len(records),'plan_sha256':b.sha256(A/'plan.json'),'scores_sha256':b.sha256(A/f'{tag}_scores.csv'),'panels_sha256':b.sha256(A/f'{tag}_panels.csv')});print('V28',tag,'complete',flush=True)

if __name__=='__main__':
 from v24_runtime import require_science
 require_science();p=argparse.ArgumentParser();p.add_argument('stage',choices=['null','treatment']);p.add_argument('--family',default='H1',choices=['H1','HCT116','HEK293T','CD4']);args=p.parse_args()
 with threadpool_limits(limits=2):
  evaluate(args.stage,args.family)
