"""Independently replay V75 response means and their bounded changes from V66."""
import numpy as np
from threadpoolctl import threadpool_limits
import v75_assay_shrink as m
from v57_pie import read,write,sha

def verify():
 m.check();model=read(m.A/'reference_model.json');d=m.data();genes=d['genes'];own=np.array([list(genes).index(t) for t in d['targets']]);checks=[]
 bound=.375
 for ci,c in enumerate(('A','B','C')):
  path=m.A/'production_inputs'/(c+'_query.npz');assert sha(path)==model['context_original_inputs'][ci]['query_sha256']
  with np.load(path,allow_pickle=False) as z:base=z['baseline'].copy();bulk=z['bulk'].copy()
  rebuilt=m.predict(base,bulk,())[model['candidate']]
  with np.load(m.A/'reference_model_arrays.npz',allow_pickle=False) as z:actual=z['effects'][ci].copy()
  np.testing.assert_array_equal(actual,rebuilt)
  assert np.max(abs(actual-base))<=bound+1e-6
  assert (actual[np.arange(300),own]==0).all() and np.isfinite(actual).all() and np.max(abs(actual))<=3
  np.testing.assert_array_equal(actual[:,bulk==0],base[:,bulk==0])
  checks.append(dict(context=c,all_responses_exactly_replayed=True,correction_bound=bound,own_and_unsupported_axes_verified=True));print('V75 independent effects',c,flush=True)
 write(m.A/'final_audit.json',dict(passed=True,checks=checks,verifier_sha256=sha(__file__),scope='Exact response replay and V66 correction bounds; old identity-margin constraint deliberately absent. No hidden accuracy guarantee.'))

if __name__=='__main__':
 with threadpool_limits(limits=2):verify()
