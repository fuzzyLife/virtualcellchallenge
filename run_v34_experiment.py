"""Sequential V34 trainable-mean flow and fixed public count diagnostics."""
import argparse,json,gc
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v34_methods as f
import v32_methods as m
import v32_runtime as runtime
from v6_methods import moments
from v25_methods import emit as old_emit
from evaluate_atlas_transfer import load_atlas_model
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from v8_methods import panel_metrics

def effect_from_values(values,query,axis,base,own,alpha):
 learned=base.copy();learned[axis]=np.clip(np.log2((values.mean(0)+.1)/(np.nan_to_num(query[axis])+.1)),-3,3);learned[own]=0;out=np.clip(base+alpha*np.clip(learned-base,-1,1),-3,3);out[own]=0;return out

def mean_screen():
 plan=f.check();runtime.install();assert not (f.A/'mean_results.json').exists();data=m.old.read_cache();genes=data['genes'];targets=data['targets'];records=[]
 for ctx in plan['contexts']:
  key=f.family(ctx['family']);net,a,rec=f.load_model(key);x,meta,measured,source=f.raw(ctx['family'],genes);query=f.normmean(x[ctx['left']]);query[~measured]=np.nan;ref=f.normmean(x[ctx['right']]);base=runtime.baseline(data,query,ctx['excluded']);own={str(t):int(np.flatnonzero(genes==t)[0]) for t in targets}
  for cell in ctx['cells']:
   t=cell['target'];ti=int(np.flatnonzero(targets==t)[0]);o=own[t];v=f.generated(net,a,x[cell['template']],query,t);sh=f.generated(net,a,x[cell['template']],query,t,True);effects={'v24':base[ti],'flow_mean25':effect_from_values(v,query,a['axis'],base[ti],o,.25),'flow_mean100':effect_from_values(v,query,a['axis'],base[ti],o,1.),'target_shuffle25':effect_from_values(sh,query,a['axis'],base[ti],o,.25)};actual=f.normmean(x[cell['actual']]);valid=measured&(ref>5);valid[o]=False;truth=np.clip(np.log2((actual+.1)/(ref+.1)),-3,3);scores={}
   for arm,e in effects.items():
    p=np.nan_to_num(query)*np.exp2(e);p[o]*=.1;p*=1e6/p.sum();lfc=np.clip(np.log2((p+.1)/(ref+.1)),-3,3);delta=lfc[valid]-truth[valid];scores[arm]={'square':float(delta@delta),'absolute':float(abs(delta).sum()),'cpm':float(np.sum((p[valid]-actual[valid])**2))}
   for arm,s in scores.items():records.append({'family':ctx['family'],'study':key,'target':t,'seed':cell['seed'],'arm':arm,'mse_ratio':s['square']/scores['v24']['square'],'mae_ratio':s['absolute']/scores['v24']['absolute'],'cpm_mse_ratio':s['cpm']/scores['v24']['cpm']})
  print('V34 mean evaluation',ctx['family'],'cases',len(ctx['cells']),flush=True);del x;gc.collect()
 study={arm:{key:{k:float(np.mean([r[k] for r in records if r['arm']==arm and r['study']==key])) for k in ['mse_ratio','mae_ratio','cpm_mse_ratio']} for key in m.OUTERS} for arm in plan['mean_arms']};summary={arm:{k:float(np.mean([r[k] for r in study[arm].values()])) for k in ['mse_ratio','mae_ratio','cpm_mse_ratio']} for arm in study};s=summary['flow_mean25'];checks={'MSE_gain':s['mse_ratio']<=.99,'MAE_gain':s['mae_ratio']<=.995,'CPM_noninferiority':s['cpm_mse_ratio']<=1.02,'worst_family':max(r['mse_ratio'] for r in study['flow_mean25'].values())<=1.02,'beats_target_shuffle':all(s[k]<summary['target_shuffle25'][k] for k in ['mse_ratio','mae_ratio'])};m.write(f.A/'mean_results.json',{'plan_sha256':m.sha(f.A/'plan.json'),'records':records,'study_summary':study,'summary':summary,'checks':checks,'passes':all(checks.values()),'scope':'Previously spent public mean-composition cases; own excluded; full flow alpha1 descriptive only. No private anchored scores.'});print(json.dumps({'summary':summary,'checks':checks},indent=2),flush=True)

def rounded(q,libs):
 q=q/q.sum(1,keepdims=True);expected=q*libs[:,None];out=np.floor(expected).astype(np.int32);frac=expected-out
 for i in range(len(q)):
  n=int(libs[i]-out[i].sum());assert n>=0
  if n:out[i,np.argpartition(frac[i],-n)[-n:]]+=1
 np.testing.assert_array_equal(out.sum(1),libs);assert out.min()>=0;return sparse.csr_matrix(out)

def count_test(stage):
 import pandas as pd
 from run_v26_experiment import distribution
 plan=f.check();runtime.install();assert not (f.A/f'{stage}_results.json').exists();data=m.old.read_cache();genes=data['genes'];targets=data['targets'];atlas,_=load_atlas_model();records=[];panels=[]
 for ctx in plan['contexts']:
  key=f.family(ctx['family']);net,a,rec=f.load_model(key);x,meta,measured,source=f.raw(ctx['family'],genes);null=stage=='null';cells=[ctx['null']] if null else ctx['cells'][:4];ids=ctx['null']['fit_controls'] if null else ctx['left'];ref=x[ctx['null']['reference'] if null else ctx['right']];p,bulk=moments(x[ids])[1:3];query=p*1e6;query[~measured]=np.nan;base=None if null else runtime.baseline(data,query,ctx['excluded']);outputs=[]
  for cell in cells:
   t=None if null else cell['target'];ti=0 if null else int(np.flatnonzero(targets==t)[0]);own=None if null else int(np.flatnonzero(genes==t)[0]);template=x[cell['template']];actual=x[cell['actual']];libs=template.sum(1).A1.astype(np.int64);effect=np.zeros(len(genes)) if null else base[ti];v=f.generated(net,a,template,query,t);new=np.zeros(len(genes)) if null else effect_from_values(v,query,a['axis'],effect,own,.25);seed=cell['seed']+ti;reference,_=old_emit(template,effect,own,seed,(libs,p,bulk),atlas);mean_count,_=old_emit(template,new,own,seed,(libs,p,bulk),atlas);q=mean_count.toarray().astype(float)/libs[:,None];neural=q.copy();neural[:,a['axis']]=v/1e6
   if own is not None:neural[:,own]=q[:,own]
   neural/=neural.sum(1,keepdims=True);shape=rounded(.75*q+.25*neural,libs);arms={'v24':reference,'flow_mean25':mean_count,'flow_shape25':shape};truth=de(actual,ref,genes);tp,tj=profile_jackknife_chunked(actual);out={'truth':tp,'truth_jack':tj};axis=np.flatnonzero(measured&~np.isin(genes,targets)&(p*1e6>5))
   for arm,count in arms.items():
    np.testing.assert_array_equal(count.sum(1).A1,libs);assert count.dtype.kind in 'iu' and (count.data>=0).all();pp,pj=profile_jackknife_chunked(count);out[arm+'_pred']=pp;out[arm+'_jack']=pj;d=de(count,ref,genes);metric=direction_metrics(d,truth,genes,-1 if own is None else own);record={'family':ctx['family'],'study':key,'target':cell['target'],'seed':cell['seed'],'arm':arm,'actual_cells':actual.shape[0],'predicted_cells':count.shape[0],'null_DE_calls':int(((d['q']<.05)&d['eligible']).sum()),'real_null_DE_calls':int(((truth['q']<.05)&truth['eligible']).sum()),**distribution(count,ref,axis),**{'real_'+k:v for k,v in distribution(actual,ref,axis).items()},**metric};records.append(record);print('V34',stage,ctx['family'],cell['target'],arm,'FID',record['fidelity'],'DE calls',record['null_DE_calls'],flush=True)
   outputs.append(out)
  cp,cj=profile_jackknife_chunked(ref)
  for arm in arms:
   pp=np.array([v[arm+'_pred'] for v in outputs]);pj=np.array([v[arm+'_jack'] for v in outputs]);tp=np.array([v['truth'] for v in outputs]);tj=np.array([v['truth_jack'] for v in outputs]);metrics={} if null else panel_metrics(pp,tp,cp,pj,tj,cj,genes,[v['target'] for v in cells],targets);rr=[v for v in records if v['family']==ctx['family'] and v['arm']==arm];avg={k:float(np.mean([v[k] for v in rr if np.isfinite(v.get(k,np.nan))])) for k in ['fidelity','reach','jaccard','nmae'] if any(np.isfinite(v.get(k,np.nan)) for v in rr)};panels.append({'family':ctx['family'],'study':key,'arm':arm,**avg,**metrics})
  del x;gc.collect()
 # Preserve undefined scorer outputs as null; they are never favorable zeros.
 def clean(v):
  if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
  if isinstance(v,list):return [clean(x) for x in v]
  if isinstance(v,(float,np.floating)) and not np.isfinite(v):return None
  return v
 m.write(f.A/f'{stage}_results.json',clean({'plan_sha256':m.sha(f.A/'plan.json'),'records':records,'panels':panels,'scope':'Raw public count metrics, no private anchor scaling. 200-cell D1 panels are diagnostics below competition400. Fixed first4 cases/context, not a random full-panel estimate.'}))

def main(stage,outer=None):
 if stage=='freeze':f.freeze()
 elif stage=='train':
  f.check()
  for key in ([outer] if outer else m.OUTERS):f.train(key)
 elif stage=='mean':mean_screen()
 else:count_test(stage)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','train','mean','null','treatment']);p.add_argument('--outer',choices=m.OUTERS);args=p.parse_args()
 with threadpool_limits(limits=2):main(args.stage,args.outer)