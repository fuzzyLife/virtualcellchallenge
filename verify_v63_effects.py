"""Verify V63 production response replay and every template-margin constraint."""
import numpy as np
from threadpoolctl import threadpool_limits
import v63_identity_guard as m
from v62_pool_cache import install
from v57_pie import sha,read,write
install()

def verify():
 m.check();m.parent.check();model=read(m.A/'reference_model.json');d=m.parent.data();genes=d['genes'];targets=d['targets'];own=np.array([list(genes).index(t) for t in targets]);panel=np.isin(genes,targets);checks=[]
 for ci,c in enumerate(('A','B','C')):
  path=m.A/'production_inputs'/(c+'_query.npz');assert sha(path)==model['context_original_inputs'][ci]['query_sha256']
  with np.load(path,allow_pickle=False) as z:q=z['query'].copy();base=z['baseline'].copy();bulk=z['bulk'].copy();unprotected=z['unprotected'].copy()
  rebuilt_source=m.parent.predict(q,base,())['assay_center50'];np.testing.assert_array_equal(rebuilt_source,unprotected)
  rebuilt,trace=m.guard(base,rebuilt_source,bulk,own,panel)
  with np.load(m.A/'reference_model_arrays.npz',allow_pickle=False) as z:actual=z['effects'][ci].copy();np.testing.assert_array_equal(rebuilt,actual)
  control=np.log1p(50000*bulk);a=(np.log1p(50000*m.probabilities(base,bulk,own))-control)[:,~panel];b=(np.log1p(50000*m.probabilities(actual,bulk,own))-control)[:,~panel];norm=np.linalg.norm(a,axis=1);u=np.divide(a,norm[:,None],out=np.zeros_like(a),where=norm[:,None]>0)
  for i in range(300):
   competitors=u[i][None,:]-u;old=competitors@a[i];new=competitors@b[i];assert np.all(new+1e-9*max(norm[i],1e-8)>=.8*np.maximum(old,0))
  checks.append(dict(context=c,source_estimator_and_all_responses_replayed=True,all_90000_template_margin_constraints_checked=True));print('V63 independent response/margin audit',c,flush=True)
 write(m.A/'final_audit.json',dict(passed=True,checks=checks,all_production_responses_replayed=True,independent_per_target_guard_inequalities=True,method_sha256=sha(m.ROOT/'v63_identity_guard.py'),verifier_sha256=sha(__file__),scope='All source effects deterministically replayed; per-target independent geometry check verifies every original-template margin. This preserves prior geometry, not hidden truth or hidden accuracy.'))

if __name__=='__main__':
 with threadpool_limits(limits=2):verify()
