"""Outcome-known V38 mean ablation: compare the learned head with zero head.

This diagnostic cannot replace a failed frozen count finalist. It tests whether
the apparent improvement is merely subtraction of a quarter of the pooled mean.
"""
import gc,json
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v38_response_hurdle as m
from run_v35_v36 import contexts,key
from run_v37_ensemble import checked_folder,aggregate
from v6_methods import moments
def main():
 folder=m.A/'specificity';folder.mkdir(exist_ok=True);record=folder/'plan.json'
 if not record.exists():
  m.write(record,{'parent_plan_sha256':m.b.sha256(m.A/'plan.json'),'script_sha256':m.b.sha256(m.ROOT/'audit_v38_head_specificity.py'),'known_outcomes':'All V38 mean/source validation outcomes known, nulls completed and first H1 treatment cases may be underway. No acceptance relaxation.','diagnostic':'Same one-quarter pooled-response subtraction with a zero learned head. Query-expression/support/own masks and clipping unchanged. No new fitted model or release arm.'})
 assert m.read(record)['script_sha256']==m.b.sha256(m.ROOT/'audit_v38_head_specificity.py')
 d=m.data();genes=d['genes'];targets=d['targets'];rows=[]
 for ctx in contexts():
  cache,_=checked_folder(ctx['family']);path=m.A/'predictions'/(ctx['family']+'.npz')
  with np.load(path,allow_pickle=False) as z:base=z['v24'].copy();learned=z['hurdle25'].copy();query=z['query'].copy();measured=z['measured'].copy()
  ds=m.descriptor(query,(key(ctx['family']),));ref=moments(sparse.load_npz(cache/'reference.npz'))[1]*1e6;pos=np.flatnonzero(np.isfinite(query)&(query>1));zero=base.copy()
  for i,cell in enumerate(ctx['cells']):
   ti=int(np.flatnonzero(targets==cell['target'])[0]);_,n=m.features(ds,ti,pos);zero[i,pos]+=.25*np.where(n>0,np.clip(-ds['pooled'][ti,pos],-1,1),0);zero[i,genes==cell['target']]=0
   actual=moments(sparse.load_npz(cache/f'{i}_actual.npz'))[1]*1e6;mask=measured&(ref>5);mask[genes==cell['target']]=False;truth=np.clip(np.log2((actual+.1)/(ref+.1)),-3,3);scores={}
   for a,eff in [('v24',base[i]),('learned',learned[i]),('zero_head',zero[i])]:
    pp=np.nan_to_num(query)*np.exp2(eff);pp[genes==cell['target']]*=.1;pp*=1e6/pp.sum();err=np.clip(np.log2((pp+.1)/(ref+.1)),-3,3)[mask]-truth[mask];scores[a]={'mse':float(err@err),'mae':float(abs(err).sum()),'cpm':float(np.sum((pp[mask]-actual[mask])**2))}
   for a in ['learned','zero_head']:rows.append({'family':ctx['family'],'study':key(ctx['family']),'target':cell['target'],'arm':a,**{k+'_ratio':scores[a][k]/scores['v24'][k] for k in scores[a]}})
  np.savez_compressed(folder/(ctx['family']+'.npz'),zero_head=zero);print('V38 head ablation',ctx['family'],flush=True);gc.collect()
 summary={a:{k:aggregate(rows,a,k) for k in ['mse_ratio','mae_ratio','cpm_ratio']} for a in ['learned','zero_head']};m.write(folder/'results.json',{'records':rows,'summary':summary,'plan_sha256':m.b.sha256(record),'scope':'Outcome-known mechanistic mean ablation, not independent confirmation or new release selection.'});print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()

