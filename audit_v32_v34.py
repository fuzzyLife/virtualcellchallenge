"""Independent audits for sequential V32-V34 pilot artifacts and boundaries."""
import argparse,json
import numpy as np
from threadpoolctl import threadpool_limits
import v32_methods as m
from functools import lru_cache
m.old.embeddings=lru_cache(maxsize=1)(m.old.embeddings)

def audit32():
 import v32_fallback_correction as c
 from run_v30_experiment import rows_for
 import v32_runtime as runtime
 m.check();runtime.install();c.install();result=json.loads((m.A/'corrected_results.json').read_text());fits=[]
 for path in sorted((m.A/'models').glob('*.npz')):
  g=c.load(path);assert all(m.family(s) not in g['excluded'] for s in g['sources']);assert not set(g['fit_targets'])&set(g['validation_targets']);d=m.paired_table(g,g['shuffled']);keep=~d['validation'];x=d['x'][keep];w=d['weight'][keep];y=d['y'][keep]/g['scale'];own=d['own'][keep];basis=g['basis'];mask=np.ones(y.shape,bool);ii=np.flatnonzero(own>=0);mask[ii,own[ii]]=False;prediction=x@g['coef']@basis.T;gradient=x.T@(w[:,None]*(prediction-y)*mask)@basis+g['regularization']*keep.sum()/x.shape[1]*g['coef'];rhs=x.T@(w[:,None]*y*mask)@basis;res=float(np.linalg.norm(gradient)/max(np.linalg.norm(rhs),1e-12));assert res<1e-6,res;fits.append({'file':path.name,'masked_gradient_relative_residual':res,'train_rows':int(keep.sum()),'sources':g['sources'],'excluded':g['excluded']});print('V32 audit fit',path.name,'gradient',res,flush=True)
 # Recompute every corrected outer summary using independently accumulated
 # float64 errors, rather than the original summarizer's float32 dot products.
 checks=[];data=m.old.read_cache()
 for r in [v for v in result['records'] if v['stage']=='outer']:
  suffix='_shuffle' if r['arm']=='identity_shuffled' else '';g=c.load(m.A/'models'/f'{r["outer"]}{suffix}.npz')
  if r['arm']=='source_copy':g['coef']*=0
  rows=rows_for(data,r['outer']);cs=[]
  for ci in np.unique(data['row_contexts'][rows]):
   pred,base,_=c.predict(g,data['profiles'][ci],r['alpha']);outside=np.ones(len(data['genes']),bool);outside[g['axis']]=False;np.testing.assert_array_equal(pred[:,outside],base[:,outside]);rr=rows[data['row_contexts'][rows]==ci];se=be=ae=ab=0.
   for ri in rr:
    ti=int(data['row_targets'][ri]);valid=data['valid'][ri].copy();valid[data['genes']==data['targets'][ti]]=False;a=(pred[ti,valid]-data['truth'][ri,valid]).astype(float);b=(base[ti,valid]-data['truth'][ri,valid]).astype(float);se+=a@a;be+=b@b;ae+=abs(a).sum();ab+=abs(b).sum()
   cs.append([se/be,ae/ab])
  vals=np.mean(cs,0);assert abs(vals[0]-r['mse_ratio'])<2e-6 and abs(vals[1]-r['mae_ratio'])<2e-6;checks.append({'outer':r['outer'],'arm':r['arm'],'replayed_MSE_ratio':float(vals[0]),'replayed_MAE_ratio':float(vals[1])})
 m.write(m.A/'independent_audit.json',{'passed':True,'corrected_results_sha256':m.sha(m.A/'corrected_results.json'),'fits':fits,'outer_replays':checks,'auditor_sha256':m.sha(m.ROOT/'audit_v32_v34.py'),'scope':'All12 fitted ridge gradients and all9 corrected outer summaries replayed; whole-family and identity splits checked. Frozen V24 fallback outside support exact.'});print('V32 independent audit passed',len(fits),'fits',len(checks),'outer comparisons',flush=True)

def audit33():
 import run_v33_experiment as t
 t.check();p=json.loads((t.A/'plan.json').read_text());r=json.loads((t.A/'results.json').read_text());assert r['plan_sha256']==m.sha(t.A/'plan.json');inventories=[]
 for path in sorted((t.A/'predictions').glob('*.npz')):
  rec=json.loads(path.with_suffix('.json').read_text());assert m.sha(path)==rec['sha256'];assert not set(rec['fit_targets'])&set(rec['validation_targets']);assert all(m.family(s)!=rec['outer'] for s in rec['fitting_sources']);assert all(m.split(v.split(':',1)[1])!=0 for v in rec['basis_training'])
  with np.load(path,allow_pickle=False) as z:
   assert z['trainx'].shape[0]<=512 and z['basis'].shape[1]==8;assert all(np.isfinite(z[arm]).all() for arm in p['arms']);q=z['queryx'];assert q.shape[1]==z['trainx'].shape[1]
  inventories.append({'outer':rec['outer'],'fitting_rows':len(rec['fit_targets']),'sources':rec['fitting_sources'],'sha256':rec['sha256']})
 assert len(inventories)==3
 for arm,v in r['summary'].items():
  for key,val in v.items():assert abs(val-np.mean([x[key] for x in r['records'] if x['arm']==arm and x['alpha']==.5]))<1e-12
 m.write(t.A/'independent_audit.json',{'passed':True,'results_sha256':m.sha(t.A/'results.json'),'models':inventories,'checks':'Checkpoint LFS SHA; all3 fitting/identity/family boundaries; program/data sizes; finite outputs and primary score aggregation. Does not independently reproduce proprietary synthetic pretraining or every GPU forward.','auditor_sha256':m.sha(m.ROOT/'audit_v32_v34.py')});print('V33 independent audit passed',flush=True)

def audit34():
 import v34_methods as f
 f.check();p=json.loads((f.A/'plan.json').read_text());records=[]
 for outer in m.OUTERS:
  a,r=f.load_training(outer);net,model,rec=f.load_model(outer);assert all(f.family(s)!=outer for s in rec['fitting_sources']);assert not set(rec['fitting_targets'])&set(rec['validation_targets']);assert all(m.split(t)!=0 for t in rec['fitting_targets'] if t!='Non-Targeting');assert all(m.split(t)==0 for t in rec['validation_targets']);assert len(rec['history'])==8;assert rec['selected_epoch']==min(rec['history'],key=lambda v:v['held_target_inference_mean_error_ratio_vs_controls'])['epoch'];assert a['axis'].shape==model['axis'].shape
  for v in r['source_records']:
   assert not set(v['control_rows'])&set(v['treatment_rows']);assert v['family']!=outer
  records.append({'outer':outer,'training_groups':len(r['groups']),'sources':rec['fitting_sources'],'selected_epoch':rec['selected_epoch'],'model_sha256':rec['sha256']})
 for name in ['mean','null','treatment']:
  r=json.loads((f.A/f'{name}_results.json').read_text());assert r['plan_sha256']==m.sha(f.A/'plan.json')
  if name=='mean':
   for arm,v in r['summary'].items():
    for k,val in v.items():assert abs(val-np.mean([s[k] for s in r['study_summary'][arm].values()]))<1e-12
  else:
   assert len(r['records'])==(15 if name=='null' else 60)
 for ctx in p['contexts']:
  n=ctx['null'];sets=[set(n[k]) for k in ['actual','reference','template']];assert all(not sets[i]&sets[j] for i in range(3) for j in range(i));assert set(n['template'])<=set(n['fit_controls']);assert not set(n['actual'])&set(n['fit_controls']);assert not set(n['reference'])&set(n['fit_controls'])
 m.write(f.A/'independent_audit.json',{'passed':True,'models':records,'results_sha256':{n:m.sha(f.A/f'{n}_results.json') for n in ['mean','null','treatment']},'checks':'All3 whole-family/target boundaries; selected epochs from source-only generated validation; control/treatment row separation; null400-cell disjointness; strict result arithmetic. Unit fixtures check own-label coupling invariance and exact generated libraries. No independent retraining or raw host-anchor reproduction.','auditor_sha256':m.sha(m.ROOT/'audit_v32_v34.py')});print('V34 independent audit passed',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('version',choices=['32','33','34']);a=p.parse_args()
 with threadpool_limits(limits=2):{'32':audit32,'33':audit33,'34':audit34}[a.version]()