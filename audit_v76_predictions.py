"""Audit V76 fold inputs, full response encoding and held-out prediction bounds."""
import numpy as np
from threadpoolctl import threadpool_limits
import v76_response_dilution as m
from v57_pie import read,write,sha

def main():
 m.check();d=m.data();genes=d['genes'];own=np.array([list(genes).index(t) for t in d['targets']]);checks=[]
 for ctx in m.contexts():
  name=ctx['family'];p=m.A/'predictions'/(name+'.npz');receipt=read(p.with_suffix('.json'));assert sha(p)==receipt['sha256']
  parent=m.ROOT/'experiments/v66_assay_strength/analysis/predictions'/(name+'.npz');assert sha(parent)==receipt['parent_sha256'];assert read(parent.with_suffix('.json'))['excluded_study']==m.study(name)
  with np.load(p,allow_pickle=False) as z:
   base=z['v66'].copy();bulk=z['bulk'].copy();p0=m.parent_prob(base,bulk,own);closed=np.clip(np.log2(np.maximum(np.divide(p0,bulk[None,:],out=np.ones_like(p0),where=bulk[None,:]>0),1e-30)),-3,3)
   genuine=m.predict(base,bulk,(m.study(name),));shuffled=m.predict(base,bulk,(m.study(name),),True)
   for arm,values in list(genuine.items())+[(a+'_shuffled',e) for a,e in shuffled.items()]:
    actual=z[arm];np.testing.assert_array_equal(actual,values);delta=abs(actual-closed);delta[np.arange(300),own]=0
    assert delta.max()<=.375+1e-6 and np.isfinite(actual).all() and np.max(abs(actual))<=3
    assert (actual[:,bulk==0]==0).all()
  checks.append(dict(source=name,inherited_whole_study_exclusion=m.study(name),all_arm_effects_exactly_replayed=True,prediction_sha256=sha(p)))
  print('V76 prediction/boundary audit',name,flush=True)
 write(m.A/'boundary_audit.json',dict(passed=True,checks=checks,auditor_sha256=sha(__file__),candidate_effects_include_own=True,candidate_emission_own_override=None))

if __name__=='__main__':
 with threadpool_limits(limits=2):main()
