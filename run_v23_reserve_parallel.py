"""Once-only raw-CD4 reserve, isolated case processes, frozen V23 numerics."""
from pathlib import Path
import json,gc,hashlib
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy.stats import wasserstein_distance
from threadpoolctl import threadpool_limits
from v23_methods import ROOT,A,NAMES,load_plan,load_model,raw,query_states,emit,logcpm,write,digest
from fit_v7_response import read_cache
from v16_methods import source_arrays
from v6_methods import moments
from v8_methods import gated,panel_metrics
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from evaluate_atlas_transfer import load_atlas_model
import build_submission as b

def one(index):
 with threadpool_limits(limits=2):
  plan=load_plan();ctx=json.loads((A/'reserve_plan.json').read_text());cell=ctx['cells'][index];data=read_cache();genes=data['genes'];targets=data['targets'];model=load_model('CD4');x,meta,measured,source=raw('CD4',genes);population=x[ctx['left']];_,p,bulk,_=moments(population);ref=x[ctx['right']];bank=source_arrays(data,['CD4']);target=cell['target'];ti=int(np.flatnonzero(targets==target)[0]);own=int(bank['own'][ti]);effect=gated(bank['base'][ti].copy(),p*1e6,own).astype(np.float32);states=query_states(population,model);lookup={v:i for i,v in enumerate(ctx['left'])};z=states[[lookup[v] for v in cell['template']]];template=x[cell['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64);actual=x[cell['actual']];assert len(set(cell['actual']))==len(set(cell['template']))==400;labels=meta.gene_target.to_numpy(str);assert (labels[cell['actual']]==target).all() and (labels[cell['template']]=='Non-Targeting').all();assert not set(cell['actual'])&(set(ctx['left'])|set(ctx['right']));truth=de(actual,ref,genes);tp,tj=profile_jackknife_chunked(actual);arrays={'truth':tp,'truth_jack':tj};atlas,_=load_atlas_model();records=[];axis=np.flatnonzero(measured&~np.isin(genes,targets)&(p*1e6>5));axis=axis[np.argsort(abs(effect[axis]),kind='stable')[-256:]];direction=effect[axis]/max(np.linalg.norm(effect[axis]),1e-9);projection=np.asarray(logcpm(actual)[:,axis]@direction).ravel();projection-=projection.mean();folder=A/'confirmation_cache';folder.mkdir(exist_ok=True);file=folder/f'CD4_{target}_{cell["seed"]}.npz';side=file.with_suffix('.json');assert not side.exists(),'Reserve already evaluated: '+target
  for name in NAMES:
   count,diag=emit(template,effect,own,cell['seed']+ti,name,(libs,p,bulk),atlas,model,ti,z);np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),libs);pp,pj=profile_jackknife_chunked(count);arrays[name+'_pred']=pp;arrays[name+'_jack']=pj;metrics=direction_metrics(de(count,ref,genes),truth,genes,own);score=np.asarray(logcpm(count)[:,axis]@direction).ravel();score-=score.mean();records.append({'family':'CD4','target':target,'seed':cell['seed'],'candidate':name,'actual_cells':400,'prediction_cells':400,'reference_cells':ref.shape[0],'count_sha256':digest(count),'shape_wasserstein':float(wasserstein_distance(projection,score)),**metrics,**diag})
  np.savez_compressed(file,**arrays);write(side,{'inputs':{'plan':b.sha256(A/'plan.json'),'names':NAMES,'model':b.sha256(A/'models/CD4.npz'),'reserve':b.sha256(A/'reserve_plan.json')},'sha256':b.sha256(file),'records':records,'original_input_record':source,'scheduler_sha256':b.sha256(Path(__file__))});print('V23 reserve case complete',target,flush=True);return target

def main():
 selection=json.loads((A/'selection.json').read_text());assert selection['new_candidate_passes'];assert not (A/'confirmation_evaluation.json').exists();ctx=json.loads((A/'reserve_plan.json').read_text());load_plan();start=A/'reserve_execution_started.json'
 if start.exists():raise FileExistsError('Reserve already opened; inspect existing caches, do not resubmit automatically')
 write(start,{'reserve_plan_sha256':b.sha256(A/'reserve_plan.json'),'selection_sha256':b.sha256(A/'selection.json'),'script_sha256':b.sha256(Path(__file__)),'policy':'Only locked mean_preserved finalist; six fixed comparison arms, no retuning.'})
 with ProcessPoolExecutor(max_workers=6) as pool:
  futures=[pool.submit(one,i) for i in range(len(ctx['cells']))]
  for future in as_completed(futures):future.result()
 data=read_cache();genes=data['genes'];targets=data['targets'];x,_,_,_=raw('CD4',genes);ref=x[ctx['right']];cp,cj=profile_jackknife_chunked(ref);scores=[];arrays=[]
 for cell in ctx['cells']:
  file=A/'confirmation_cache'/f'CD4_{cell["target"]}_{cell["seed"]}.npz';rec=json.loads(file.with_suffix('.json').read_text());assert b.sha256(file)==rec['sha256'];scores.extend(rec['records'])
  with np.load(file,allow_pickle=False) as z:arrays.append({k:z[k].copy() for k in z.files})
 tp=np.array([v['truth'] for v in arrays]);tj=np.array([v['truth_jack'] for v in arrays]);panels=[]
 for name in NAMES:
  pp=np.array([v[name+'_pred'] for v in arrays]);pj=np.array([v[name+'_jack'] for v in arrays]);m=panel_metrics(pp,tp,cp,pj,tj,cj,genes,[c['target'] for c in ctx['cells']],targets);f=pd.DataFrame([r for r in scores if r['candidate']==name]);panels.append({'family':'CD4','seed':ctx['cells'][0]['seed'],'candidate':name,**f.mean(numeric_only=True).to_dict(),**m})
 pd.DataFrame(scores).to_csv(A/'confirmation_scores.csv',index=False);f=pd.DataFrame(panels);f.to_csv(A/'confirmation_panels.csv',index=False);fields=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio','shape_wasserstein'];summary=f.groupby('candidate')[fields].mean();summary.to_csv(A/'confirmation_summary.csv');write(A/'confirmation_evaluation.json',{'complete':True,'populations':len(scores),'case_count':len(ctx['cells']),'plan_sha256':b.sha256(A/'plan.json'),'scores_sha256':b.sha256(A/'confirmation_scores.csv'),'panels_sha256':b.sha256(A/'confirmation_panels.csv'),'reserve_plan_sha256':b.sha256(A/'reserve_plan.json'),'scheduler_sha256':b.sha256(Path(__file__)),'scope':'New raw-cell donor reserve from an existing study; no private anchor/score estimate.'});print(summary.to_string(),flush=True)
if __name__=='__main__':main()
