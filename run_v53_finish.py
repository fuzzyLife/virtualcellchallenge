"""Verify linked V52 labels, then evaluate and audit completed V53 models."""
from pathlib import Path
import subprocess,sys,time
import numpy as np
import run_v53_joint_guard as guard
import v52_joint_response as m
ROOT=m.ROOT;A=guard.A;BASE=guard.BASE

def run(*args):
 print('V53 finish:',*args,flush=True);subprocess.run([sys.executable,*args],cwd=ROOT,check=True)

def linked_labels():
 p=m.read(BASE/'label_audit.json');assert p['passed'];links={}
 for source in m.SOURCES:
  old=m.read(BASE/'labels'/(source+'.json'));new=m.read(A/'labels'/(source+'.json'));linked=dict(new)
  linked.pop('linked_npz');linked.pop('parent_label_plan_sha256');linked['plan_sha256']=old['plan_sha256'];assert linked==old
  path=BASE/'labels'/(source+'.npz');assert m.b.sha256(path)==old['sha256'];links[source]=old['sha256']
 m.write(A/'label_audit.json',dict(passed=True,parent_raw_cell_audit_sha256=m.b.sha256(BASE/'label_audit.json'),linked_label_sha256=links,
  scope='Reuses independent V52 raw-cell pairwise replay after verifying identical label bytes and cell/view metadata. No labels recomputed or modified.',plan_sha256=m.b.sha256(A/'plan.json'),auditor_sha256=m.b.sha256(Path(__file__)),D3_opened=False))

def feature_contract():
 guard.install();proof=[]
 for source in ['H1','HCT116','CD4_Rest_D2']:
  outer=m.GROUP[source];z=m.labels(source);ti=int(np.flatnonzero(z['valid'].any(1))[0]);pos=np.flatnonzero(z['valid'][ti]);x,_=m.features(source,outer,ti,pos)
  with np.load(A/'models'/outer/'normalizer.npz',allow_pickle=False) as norm:med=norm['median'];scale=norm['scale']
  active=np.isfinite(scale);out=guard.transform(x,med,scale);n=x.shape[1]
  assert np.isfinite(out).all();np.testing.assert_array_equal(out[:,np.r_[~active,~active]],0)
  proof.append(dict(outer=outer,source=source,disabled_numeric_channels=int((~active).sum()),raw_numeric_channels=n,
   disabled_channels_observed_in_query=int(((~active)&np.isfinite(x).any(0)).sum()),disabled_channels_stay_zero=True))
 m.write(A/'feature_contract_audit.json',dict(passed=True,cases=proof,plan_sha256=m.b.sha256(A/'plan.json'),auditor_sha256=m.b.sha256(Path(__file__))))

def main():
 start=time.monotonic()
 while not all((A/'models'/outer/'complete.json').exists() for outer in ['H1','Orion','CD4','all']):
  if time.monotonic()-start>2400:raise TimeoutError('V53 training not complete; inspect train.log')
  time.sleep(10)
 linked_labels();run('run_v53_joint_guard.py','evaluate');run('run_v53_joint_guard.py','audit','--kind','boundaries');feature_contract();run('run_v53_joint_guard.py','audit','--kind','artifacts');run('document_v53_joint_guard.py')
 print('V53 completed and audited',flush=True)
if __name__=='__main__':main()
