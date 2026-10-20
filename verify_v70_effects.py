"""Replay V70 effects and independently check all V66 identity-margin bounds."""
import numpy as np
from threadpoolctl import threadpool_limits
import v70_count_teacher as m
from v63_identity_guard import probabilities
from v57_pie import read,write,sha

def verify():
 m.check();model=read(m.A/'reference_model.json');d=m.data();genes=d['genes'];own=np.array([list(genes).index(t) for t in d['targets']]);panel=np.isin(genes,d['targets']);checks=[]
 for ci,c in enumerate(('A','B','C')):
  path=m.A/'production_inputs'/(c+'_query.npz');assert sha(path)==model['context_original_inputs'][ci]['query_sha256']
  with np.load(path,allow_pickle=False) as z:base=z['baseline'].copy();bulk=z['bulk'].copy()
  rebuilt=m.predict(base,bulk,())[model['candidate']]
  with np.load(m.A/'reference_model_arrays.npz',allow_pickle=False) as z:actual=z['effects'][ci].copy()
  np.testing.assert_array_equal(actual,rebuilt)
  control=np.log1p(50000*bulk);aa=(np.log1p(50000*probabilities(base,bulk,own))-control)[:,~panel];bb=(np.log1p(50000*probabilities(actual,bulk,own))-control)[:,~panel]
  norm=np.linalg.norm(aa,axis=1);unit=np.divide(aa,norm[:,None],out=np.zeros_like(aa),where=norm[:,None]>0)
  for i in range(300):
   axes=unit[i][None,:]-unit;old=axes@aa[i];new=axes@bb[i]
   assert np.all(new+1e-9*max(norm[i],1e-8)>=.8*np.maximum(old,0)),(c,i)
  checks.append(dict(context=c,all_responses_exactly_replayed=True,all_90000_margin_constraints_verified=True));print('V70 independent effects',c,flush=True)
 write(m.A/'final_audit.json',dict(passed=True,checks=checks,verifier_sha256=sha(__file__),scope='Exact response replay and V66 geometry constraints; not a hidden accuracy claim.'))

if __name__=='__main__':
 with threadpool_limits(limits=2):verify()
