"""Audit V62 source-template formula, study exclusion and all production effects."""
import gc
import numpy as np
from v57_pie import read,write,sha
import v62_conserved_response as m
from v62_pool_cache import install
install()

def verify():
 model=read(m.A/'reference_model.json');plan=read(m.A/'production_plan.json');candidate=plan['candidate'];checks=[]
 bank=m.bank(());mask=bank['ok'] & ((bank['ok'].sum(2)>0).sum(0)>=2)[None,:,None]
 # Independent leave-one-target-out reconstruction, different from broadcast total subtraction.
 for fi in range(len(m.FAMILIES)):
  for ti in (0,71,149,299):
   eligible=np.flatnonzero(np.arange(300)!=ti);z=mask[fi,eligible];n=z.sum(0);t=np.divide(np.sum(np.where(z,bank['y'][fi,eligible],0),axis=0),n,out=np.zeros(len(n)),where=n>=20)
   np.testing.assert_allclose(t,bank['source_templates'][fi,ti],atol=2e-6,rtol=2e-5)
 for ci,c in enumerate(('A','B','C')):
  path=m.A/'production_inputs'/(c+'_query.npz');assert sha(path)==model['context_original_inputs'][ci]['query_sha256']
  with np.load(path,allow_pickle=False) as z:q=z['query'].copy();base=z['baseline'].copy()
  rebuilt=m.predict(q,base,())[candidate]
  with np.load(m.A/'reference_model_arrays.npz',allow_pickle=False) as z:np.testing.assert_array_equal(rebuilt,z['effects'][ci])
  checks.append(dict(context=c,every_response_replayed=True,coordinates=int(rebuilt.size)));print('V62 response replay',c,flush=True)
 write(m.A/'final_audit.json',dict(passed=True,independent_leave_one_target_out_formula_checks=20,all_production_response_replay=True,checks=checks,method_sha256=sha(m.ROOT/'v62_conserved_response.py'),verifier_sha256=sha(__file__),scope='Independent template arithmetic plus deterministic all-response replay; reuses the frozen pooling kernel and does not establish hidden accuracy.'))

if __name__=='__main__':
 from threadpoolctl import threadpool_limits
 with threadpool_limits(limits=2):verify()
