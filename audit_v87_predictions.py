"""Replay V87 predicted profile margins and raw public metric aggregation independently."""
from pathlib import Path
import argparse
import numpy as np
from threadpoolctl import threadpool_limits
import v87_conditional_adapter as m
from v57_pie import read,write,sha
from run_v35_v36 import contexts,key

def predictions():
 m.check();d=m.data();own=np.array([list(d['genes']).index(t) for t in d['targets']]);axis=~np.isin(d['genes'],d['targets']);records=[]
 for ctx in contexts():
  name=ctx['family'];p=m.A/'predictions'/(name+'.npz');assert sha(p)==read(p.with_suffix('.json'))['sha256']
  with np.load(p,allow_pickle=False) as z:a={k:z[k].copy() for k in z.files}
  source=m.ROOT/'experiments/v66_assay_strength/analysis/predictions'/(name+'.npz')
  with np.load(source,allow_pickle=False) as z:
   for field,old in [('v66','assay75'),('query','query'),('bulk','bulk')]:np.testing.assert_array_equal(a[field],z[old])
  def delta(effect):
   p=np.exp2(np.asarray(effect,float))*a['bulk'];p[np.arange(300),own]=.1*a['bulk'][own];p/=p.sum(1,keepdims=True);return (np.log1p(50000*p)-np.log1p(50000*a['bulk']))[:,axis]
  base=delta(a['v66']);norm=np.linalg.norm(base,axis=1);unit=base/np.maximum(norm[:,None],1e-30);s=base@unit.T;floor=.8*np.maximum(np.diag(s)[:,None]-s,0);np.fill_diagonal(floor,-np.inf)
  for arm in m.ARMS:
   assert a[arm].shape==(300,len(d['genes'])) and np.isfinite(a[arm]).all();assert (a[arm]>=-3).all() and (a[arm]<=3).all();assert (a[arm][np.arange(300),own]==0).all();pred=delta(a[arm]);ss=pred@unit.T;observed=np.diag(ss)[:,None]-ss;assert np.all(observed+1e-9*np.maximum(norm,1e-8)[:,None]>=floor)
   tr=read(p.with_suffix('.json'))['trace'][arm]
   if 'routing' in tr:assert all(m.GROUP[n]!=key(name) for n in tr['routing']['source_names'])
   records.append(dict(context=name,arm=arm,all_90000_margin_constraints_passed=True,changed_targets=int(np.any(a[arm]!=a['v66'],axis=1).sum())))
  print('V87 independent prediction audit',name,flush=True)
 write(m.A/'prediction_audit.json',dict(passed=True,records=records,own_knockdown_applied_once_in_profile_replay=True,V66_inputs_exact=True,script_sha256=sha(__file__)))

def metrics():
 from evaluate_v77 import profiles
 from audit_transfer_validation import pds
 result=read(m.A/'mean_results.json');d=m.data();own=np.array([list(d['genes']).index(t) for t in d['targets']]);records=[]
 for r in result['records']:
  p=m.A/'predictions'/(r['source']+'.npz');tp=m.ROOT/'experiments/v77_common_response/analysis/truth_bank'/(r['source']+'.npz');assert sha(tp)==read(tp.with_suffix('.json'))['sha256']
  with np.load(tp,allow_pickle=False) as z:b={k:z[k].copy() for k in z.files}
  with np.load(p,allow_pickle=False) as z:
   for arm in ('v66',)+m.ARMS:
    f=profiles(z[arm],z['bulk'],own);e=[]
    for i,ti in enumerate(b['target_indices']):
     mask=b['measured'].copy();mask[own[ti]]=False;e.append(np.mean(np.square(f[ti,mask]-b['truth'][i,mask])))
    np.testing.assert_allclose(e,r['profile_errors'][arm],rtol=1e-12,atol=1e-15);np.testing.assert_allclose(np.mean(e),r['scores'][arm]['profile_MSE'],rtol=1e-12)
    # Independent explicit ranking loop over cosine distances, with half credit for ties.
    axis=b['measured']&~np.isin(d['genes'],d['targets']);x=np.asarray(f[b['target_indices']][:,axis],float)-b['control'][axis];y=np.asarray(b['truth'][:,axis],float)-b['control'][axis];den=np.linalg.norm(x,axis=1)[:,None]*np.linalg.norm(y,axis=1)[None,:];dist=1-np.divide(x@y.T,den,out=np.zeros((len(x),len(y))),where=den>0)
    rank=[]
    for i,row in enumerate(dist):
     less=int(np.sum(row<row[i]));equal=int(np.sum(row==row[i]))-1;rank.append(1-(less+.5*equal)/(len(row)-1))
    got=float(np.mean(rank));np.testing.assert_allclose(got,r['scores'][arm]['pds'],rtol=0,atol=1e-12)
  records.append(dict(context=r['source'],pairs=r['pairs'],profile_errors_replayed=True,pds_replayed=True))
 for arm in m.ARMS:
  vals=[]
  for study in ['H1','Orion','CD4']:
   rows=[r for r in result['records'] if r['study']==study];v=np.mean([r['scores'][arm]['profile_MSE']/r['scores']['v66']['profile_MSE'] for r in rows]);np.testing.assert_allclose(v,result['studies'][arm][study]['profile_ratio'],rtol=1e-12);vals.append(v)
  np.testing.assert_allclose(np.mean(vals),result['summary'][arm]['profile_ratio'],rtol=1e-12)
 s=result['summary']['conditional25'];checks=dict(profile=s['profile_ratio']<=.995,worst_study=max(r['profile_ratio'] for r in result['studies']['conditional25'].values())<=1.03,pds=s['pds_change']>=-.005,beats_constant=s['profile_ratio']<result['summary']['constant25']['profile_ratio'],beats_shuffle=s['profile_ratio']<result['summary']['conditional25_shuffled']['profile_ratio']);assert checks==result['checks'] and all(checks.values())==result['passes']
 write(m.A/'metric_audit.json',dict(passed=True,records=records,equal_study_aggregation_replayed=True,qualification_replayed=True,scope='Profile errors independently recomputed; PDS independently replayed with explicit tied-distance ranks. No anchored-score equivalence claim.',script_sha256=sha(__file__)))
 print('V87 independent metric aggregation audit passed',flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['predictions','metrics']);a=p.parse_args()
 with threadpool_limits(limits=2):{'predictions':predictions,'metrics':metrics}[a.stage]()
