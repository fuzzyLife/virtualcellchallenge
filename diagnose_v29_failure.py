"""Post-rejection H1 mechanism diagnostics; no new candidate or score selection.

Use already-spent public treatment cases to decompose condition-dependent
logits and describe batch/depth/variance differences. Actual metadata is used
only for this retrospective diagnostic, never in prediction or release.
"""
from pathlib import Path
import json,gc
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy import sparse
from threadpoolctl import threadpool_limits
from v29_methods import ROOT,A,load_plan,load_model,write
from v26_methods import raw,moments,baseline,state_values,descriptor,decode
from v22_methods import logcpm
from fit_v7_response import read_cache
import build_submission as b

def stats(count):
 z=logcpm(count);n=len(count.indptr)-1;mu=np.asarray(z.mean(0)).ravel();var=np.maximum(np.asarray(z.power(2).mean(0)).ravel()-mu**2,0)*n/max(n-1,1);detect=np.asarray((count>0).mean(0)).ravel();return mu,var,detect

def main():
 decision=json.loads((A/'release_decision.json').read_text());assert not decision['released'];plan=load_plan();ctx=next(c for c in plan['contexts'] if c['family']=='H1');data=read_cache();genes=data['genes'];targets=data['targets'];x,meta,measured,_=raw('H1',genes);ids=ctx['left'];p,bulk=moments(x[ids])[1:3];query=p*1e6;query[~measured]=np.nan;effects=baseline(data,query,ctx['excluded']);allctrl=np.array(ctx['left']+ctx['right']);allstats=stats(x[allctrl]);axis=measured&~np.isin(genes,targets)&(p*1e6>5);rows=[];variance=[]
 with np.load(ROOT/'experiments/v27_control_anchored/analysis/training/H1.npz',allow_pickle=False) as z:source_cond=z['conditions'].copy()
 norms=np.linalg.norm(source_cond,axis=1);norms=norms[norms>0]
 for panel in [False,True]:
  model=load_model('H1',panel);ss=state_values(x,meta,ids,model['state_axis'],model['state_basis'],model['state_scale']);lookup={int(v):i for i,v in enumerate(ids)}
  for cell in ctx['cells']:
   ti=int(np.flatnonzero(targets==cell['target'])[0]);seed=cell['seed']+ti;state=ss[[lookup[int(v)] for v in cell['template']]];condition=descriptor(effects[ti],p,model['geometry']);cc=np.repeat(condition[None,:],400,axis=0);latent=np.random.default_rng(np.random.SeedSequence([seed,27027])).normal(size=(400,16)).astype(np.float32);delta=np.clip(decode(model,latent,state,cc)-decode(model,latent,state,np.zeros_like(cc)),-1,1);mask=model['output_mask'].copy();mask[genes==cell['target']]=False;delta=delta[:,mask];mean=delta.mean(0);centered=delta-mean;total=float(np.mean(delta**2));interaction=float(np.mean(centered**2));rows.append({'target':cell['target'],'seed':cell['seed'],'model':'panel' if panel else 'masked','query_condition_norm':float(np.linalg.norm(condition)),'median_source_condition_norm':float(np.median(norms)),'condition_delta_RMS':float(np.sqrt(total)),'cell_interaction_RMS':float(np.sqrt(interaction)),'cell_interaction_energy_fraction':interaction/max(total,1e-30),'scope':'Constant per-gene condition offsets are largely removed by mean projection; centered component can modify cell variation.'})
 for cell in ctx['cells']:
  actual=x[cell['actual']];actualstats=stats(actual);batches=meta['sample'].to_numpy(str);freq=pd.Series(batches[cell['actual']]).value_counts(normalize=True);parts=[];weights=[];nctrl=[];matched_ids=[]
  for name,w in freq.items():
   use=allctrl[batches[allctrl]==name]
   if len(use)<5:continue
   parts.append(stats(x[use]));weights.append(float(w));nctrl.append(len(use));matched_ids.extend(use.tolist())
  if len(parts)!=len(freq):continue
  w=np.array(weights);w/=w.sum();mu=np.sum(w[:,None]*np.array([s[0] for s in parts]),axis=0);within=np.sum(w[:,None]*np.array([s[1] for s in parts]),axis=0);between=np.sum(w[:,None]*(np.array([s[0] for s in parts])-mu)**2,axis=0);vv=within+between;dd=np.sum(w[:,None]*np.array([s[2] for s in parts]),axis=0);libs=actual.sum(1).A1;variance.append({'target':cell['target'],'seed':cell['seed'],'treated_batches':len(freq),'pooled_NTC_batches':len(set(batches[allctrl])),'matched_control_cells':sum(nctrl),'median_treated_depth':float(np.median(libs)),'median_pooled_control_depth':float(np.median(x[allctrl].sum(1).A1)),'median_matched_control_depth':float(np.median(x[matched_ids].sum(1).A1)),'median_log_variance_ratio_actual_pooled_NTC':float(np.median(np.log((actualstats[1][axis]+.01)/(allstats[1][axis]+.01)))),'median_log_variance_ratio_actual_matched_NTC':float(np.median(np.log((actualstats[1][axis]+.01)/(vv[axis]+.01)))),'mean_between_batch_fraction_matched_NTC':float(np.mean(between[axis]/np.maximum(vv[axis],1e-12))),'detection_MAE_actual_pooled_NTC':float(abs(actualstats[2][axis]-allstats[2][axis]).mean()),'detection_MAE_actual_matched_NTC':float(abs(actualstats[2][axis]-dd[axis]).mean())})
 pd.DataFrame(rows).to_csv(A/'condition_mechanism_diagnostic.csv',index=False);pd.DataFrame(variance).to_csv(A/'batch_variance_diagnostic.csv',index=False);write(A/'failure_diagnostic.json',{'release_decision_sha256':b.sha256(A/'release_decision.json'),'plan_sha256':b.sha256(A/'plan.json'),'condition_rows_sha256':b.sha256(A/'condition_mechanism_diagnostic.csv'),'batch_rows_sha256':b.sha256(A/'batch_variance_diagnostic.csv'),'diagnostician_sha256':b.sha256(Path(__file__)),'source_condition_norm_median':float(np.median(norms)),'scope':'Descriptive post-rejection analysis of spent H1 biology. Batch-matched actual metadata is not available in query inference and is not used in a candidate. No reserve opening, new eligible model, threshold change or score reselection.'});print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(variance).to_string(index=False));print('V29 descriptive failure diagnostics complete',flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
