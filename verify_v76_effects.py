"""Replay all V76 effects and verify composition/dose/step bounds."""
import numpy as np
from threadpoolctl import threadpool_limits
import v76_response_dilution as m
from v57_pie import read,write,sha

def verify():
 m.check();model=read(m.A/'reference_model.json');d=m.data();genes=d['genes'];own=np.array([list(genes).index(t) for t in d['targets']]);checks=[]
 parent=m.ROOT/'experiments/v66_assay_strength/analysis/reference_model_arrays.npz';assert sha(parent)==model['parent_V66_arrays_sha256']
 with np.load(parent,allow_pickle=False) as z:parent_effects=z['effects'].copy()
 for ci,c in enumerate(('A','B','C')):
  p=m.A/'production_inputs'/(c+'_query.npz');assert sha(p)==model['context_original_inputs'][ci]['query_sha256']
  with np.load(p,allow_pickle=False) as z:base=z['baseline'].copy();bulk=z['bulk'].copy()
  np.testing.assert_array_equal(base,parent_effects[ci]);rebuilt=m.predict(base,bulk,())[model['candidate']]
  with np.load(m.A/'reference_model_arrays.npz',allow_pickle=False) as z:actual=z['effects'][ci].copy()
  np.testing.assert_array_equal(actual,rebuilt);assert np.isfinite(actual).all() and np.max(abs(actual))<=3
  p0=m.parent_prob(base,bulk,own);closed=np.clip(np.log2(np.maximum(np.divide(p0,bulk[None,:],out=np.ones_like(p0),where=bulk[None,:]>0),1e-30)),-3,3)
  change=abs(actual-closed);change[np.arange(300),own]=0;assert change.max()<=.375+1e-6
  assert (actual[:,bulk==0]==0).all()
  if model['candidate'].startswith('mixture'):assert actual[np.arange(300),own].max()<=1e-6
  checks.append(dict(context=c,all_responses_exactly_replayed=True,nonown_correction_bound=.375,full_own_effects_explicit=True,unexpressed_genes_zero=True));print('V76 independent effects',c,flush=True)
 write(m.A/'final_audit.json',dict(passed=True,checks=checks,verifier_sha256=sha(__file__),scope='Exact response replay; full encoded own effects are applied once with no legacy own override. No hidden accuracy guarantee.'))

if __name__=='__main__':
 with threadpool_limits(limits=2):verify()
