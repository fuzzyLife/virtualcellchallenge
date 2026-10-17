"""Independent label replay, source-boundary mutation and V52 artifact audit."""
import argparse,gc,json
from pathlib import Path
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
import v52_joint_response as m

def label_replay(selected=None):
 results=[]
 for source in ([selected] if selected else m.SOURCES):
  z=m.labels(source);rec=m.read(m.A/'labels'/(source+'.json'))
  assert len(rec['groups'])==int(z['valid'].any(1).sum());assert z['mean'].shape==z['rank'].shape==(300,2,len(z['genes']))
  assert np.isfinite(z['mean']).all() and np.isfinite(z['rank']).all() and abs(z['rank']).max()<=.5+1e-7
  for group in rec['groups']:
   ti=group['target_index'];own=int(np.flatnonzero(z['genes']==group['target'])[0]);assert not z['valid'][ti,own]
   left=set();right=set();ctrl0=set();ctrl1=set()
   for lane in group['lanes']:
    left.update(lane['treated_views'][0]);right.update(lane['treated_views'][1]);ctrl0.update(lane['control_views'][0]);ctrl1.update(lane['control_views'][1])
   assert not left&right and not ctrl0&ctrl1 and not (left|right)&(ctrl0|ctrl1)
  group=rec['groups'][0];ti=group['target_index'];pos=np.flatnonzero(z['valid'][ti])[:32]
  x,meta,measured,original=m.raw(source,z['genes']);labels=meta.gene_target.to_numpy(str);batches=meta['sample'].to_numpy(str)
  # Independent pair-comparison implementation, without scipy rankdata.
  for view in range(2):
   total=sum(len(l['treated_views'][view]) for l in group['lanes']);mean=np.zeros(len(pos));control=mean.copy();rank=mean.copy()
   for lane in group['lanes']:
    tt=np.array(lane['treated_views'][view]);cc=np.array(lane['control_views'][view]);assert np.all(labels[tt]==group['target']) and np.all(labels[cc]=='Non-Targeting')
    assert np.all(batches[tt]==lane['lane']) and np.all(batches[cc]==lane['lane'])
    tx=x[tt].toarray().astype(float);cx=x[cc].toarray().astype(float);tx=tx[:,pos]*(1e6/tx.sum(1))[:,None];cx=cx[:,pos]*(1e6/cx.sum(1))[:,None]
    weight=len(tt)/total;diff=tx[:,None,:]-cx[None,:,:]
    rank+=weight*((diff>0).mean((0,1))+.5*(diff==0).mean((0,1))-.5);mean+=weight*tx.mean(0);control+=weight*cx.mean(0)
   truth=np.clip(np.log2((mean+.1)/(control+.1)),-3,3)
   np.testing.assert_allclose(z['mean'][ti,view,pos],truth,atol=2e-7,rtol=1e-6)
   np.testing.assert_allclose(z['rank'][ti,view,pos],rank,atol=2e-8,rtol=1e-6)
  results.append(dict(source=source,groups=len(rec['groups']),independent_replay_target=group['target'],genes=len(pos),views=2,passed=True));del x;gc.collect()
 m.write(m.A/(('label_audit_'+selected+'.json') if selected else 'label_audit.json'),dict(passed=True,sources=results,plan_sha256=m.b.sha256(m.A/'plan.json'),auditor_sha256=m.b.sha256(Path(__file__)),D3_opened=False))
 print('V52 independent labels passed',flush=True)

def boundaries():
 source='HCT116';outer='H1';z=m.labels(source);ti=int(np.flatnonzero(z['valid'].any(1))[0]);pos=np.flatnonzero(z['valid'][ti])[:64]
 before,prior=m.features(source,outer,ti,pos)
 from predict_v52_joint_response import query_state,query_features
 replay,rprior=query_features(query_state(z['query'],outer,m.GROUP[source]),ti,pos)
 np.testing.assert_array_equal(before,replay);np.testing.assert_array_equal(prior,rprior)
 bank=m.labels('H1');original=bank['rank'].copy();d=m.data();rows=d['row_families']=='H1';saved=d['truth'][rows].copy()
 try:
  bank['rank'][:]=.49;d['truth'][rows]=2.9;m.context_features.cache_clear();m.old.bank.cache_clear()
  after,other=m.features(source,outer,ti,pos)
  np.testing.assert_array_equal(before,after);np.testing.assert_array_equal(prior,other)
 finally:
  bank['rank'][:]=original;d['truth'][rows]=saved;m.context_features.cache_clear();m.old.bank.cache_clear()
 allowed=[study for study,*_ in m.donor_rank('H1','Orion')];assert allowed==['CD4']
 m.write(m.A/'boundary_audit.json',dict(passed=True,mutation='Replace excluded H1 mean response with 2.9 and rank response with .49; Orion fitting features and priors unchanged',tested_source=source,outer=outer,generic_inference_feature_parity=True,rank_donor_studies=allowed,features=len(before),plan_sha256=m.b.sha256(m.A/'plan.json'),auditor_sha256=m.b.sha256(Path(__file__)),D3_opened=False))
 print('V52 excluded-label mutation passed',flush=True)

def artifacts():
 m.check();proof={};records=[]
 for outer in [*m.OUTERS,'all']:
  folder=m.A/'models'/outer;complete=m.read(folder/'complete.json')
  for name,h in complete['models'].items():assert m.b.sha256(folder/name)==h;proof[str((folder/name).relative_to(m.ROOT))]=h
  assert m.b.sha256(folder/'linear.joblib')==complete['linear_sha256']
  table=m.read(folder/'training_table.json');assert all(r['study']!=outer for r in table['groups']) if outer!='all' else True
  for r in table['groups']:assert r['validation']==(m.target_split(r['target'])==0)
  records.append(dict(outer=outer,groups=len(table['groups']),training_pairs=table['training_pairs'],validation_pairs=table['validation_pairs']))
 for source in m.SOURCES:
  p=m.A/'predictions'/(source+'.npz');r=m.read(p.with_suffix('.json'));assert m.b.sha256(p)==r['sha256'] and r['excluded_study']==m.GROUP[source]
  with np.load(p,allow_pickle=False) as z:
   assert z['mean'].shape==z['rank'].shape==z['valid'].shape==(300,18533)
   assert np.isfinite(z['mean']).all() and np.isfinite(z['rank']).all();assert abs(z['mean']).max()<=3 and abs(z['rank']).max()<=.5
  proof[str(p.relative_to(m.ROOT))]=r['sha256']
 from predict_v52_joint_response import predict
 z=m.labels('H1');ti=int(np.flatnonzero(z['valid'].any(1))[0]);t=z['targets'][ti]
 api=predict(z['query'],outer='H1',recipient_study='H1',targets=[t],device='cpu')
 with np.load(m.A/'predictions/H1.npz',allow_pickle=False) as saved:
  mask=z['valid'][ti]&api['valid'][0]
  np.testing.assert_allclose(api['mean'][0,mask],saved['mean'][ti,mask],atol=3e-6,rtol=1e-5)
  np.testing.assert_allclose(api['rank'][0,mask],saved['rank'][ti,mask],atol=3e-6,rtol=1e-5)
 assert not api['valid'][0,z['genes']==t].any()
 for n in ['label_audit.json','boundary_audit.json']:assert m.read(m.A/n)['passed']
 m.write(m.A/'final_audit.json',dict(passed=True,generic_inference_matches_saved_heldout_predictions=True,models=proof,training=records,plan_sha256=m.b.sha256(m.A/'plan.json'),auditor_sha256=m.b.sha256(Path(__file__)),D3_opened=False,
  scope='Model/source-label/inference integrity; does not certify private generalization or a count decoder.'))
 print('V52 final artifacts passed',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['labels','boundaries','artifacts']);p.add_argument('--source',choices=m.SOURCES);args=p.parse_args()
 with threadpool_limits(limits=2):{'labels':lambda:label_replay(args.source),'boundaries':boundaries,'artifacts':artifacts}[args.stage]()
