"""Replay all V65 production responses and independently check the prior-margin bounds."""
import numpy as np
from threadpoolctl import threadpool_limits
import v65_closed_response as m
from v63_identity_guard import probabilities
from v57_pie import read,write,sha

def verify():
 m.check();model=read(m.A/'reference_model.json');d=m.old.data();own=np.array([list(d['genes']).index(t) for t in d['targets']]);panel=np.isin(d['genes'],d['targets']);checks=[]
 for ci,c in enumerate(('A','B','C')):
  path=m.A/'production_inputs'/(c+'_query.npz');assert sha(path)==model['context_original_inputs'][ci]['query_sha256']
  with np.load(path,allow_pickle=False) as z:q=z['query'].copy();base=z['baseline'].copy();bulk=z['bulk'].copy()
  rebuilt,_=m.predict(q,base,bulk,())
  with np.load(m.A/'reference_model_arrays.npz',allow_pickle=False) as z:actual=z['effects'][ci].copy()
  np.testing.assert_array_equal(rebuilt[model['candidate']],actual)
  ctrl=np.log1p(50000*bulk);a=(np.log1p(50000*probabilities(base,bulk,own))-ctrl)[:,~panel];b=(np.log1p(50000*probabilities(actual,bulk,own))-ctrl)[:,~panel];norm=np.linalg.norm(a,axis=1);u=np.divide(a,norm[:,None],out=np.zeros_like(a),where=norm[:,None]>0)
  for i in range(300):
   competitors=u[i][None,:]-u;old=competitors@a[i];new=competitors@b[i];assert np.all(new+1e-9*max(norm[i],1e-8)>=.8*np.maximum(old,0))
  checks.append(dict(context=c,source_and_all_responses_replayed=True,all_90000_margin_constraints_verified=True));print('V65 independent response audit',c,flush=True)
 write(m.A/'final_audit.json',dict(passed=True,checks=checks,all_production_responses_replayed=True,verifier_sha256=sha(__file__),scope='Deterministic response replay and independent margin inequalities preserve V24 model geometry; no hidden accuracy guarantee.'))
if __name__=='__main__':
 with threadpool_limits(limits=2):verify()
