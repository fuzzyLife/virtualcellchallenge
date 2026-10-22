"""Numerical and provenance audit of final-panel fitted effects; no final accuracy claim."""
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
import final_response_methods as f
import v86_joint_posterior as original
from v57_pie import read,write,sha
R=f.R;A=f.A
def main():
 f.check();post=f.posterior();bank=post.bank('production');y,v,ok=post.observations(bank);active=np.flatnonzero(ok.all(1));assert len(active)>10;ids=active[np.linspace(0,len(active)-1,min(128,len(active))).astype(int)];cov=.1**2*post.kernel(bank['similarity'],.2);fast=post.component_likelihood(y[ids],v[ids],cov);slow=original.component_likelihood(y[ids],v[ids],cov)
 for x,z in zip(fast,slow):np.testing.assert_allclose(x,z,rtol=1e-10,atol=1e-10)
 priors=[]
 for shuffled in [False,True]:
  model=post.model('production',shuffled);assert (model['weights']>0).all() and np.isfinite(model['weights']).all();np.testing.assert_allclose(model['weights'].sum(),1,atol=1e-12);rec=read(A/'models'/('production_shuffled.json' if shuffled else 'production.json'));assert rec['source_cache_sha256']==sha(A/'source_cache.npz');assert rec['plan_sha256']==sha(A/'plan.json');np.testing.assert_array_equal(model['axis'],bank['axis']);priors.append(dict(shuffled=shuffled,iterations=len(rec['em_history']),last_maximum_weight_change=rec['em_history'][-1]['maximum_weight_change'],weights_sum=float(model['weights'].sum()),finite_update_budget=100,convergence_claim=False))
 r=read(A/'reference_model.json');assert sha(A/'reference_model_arrays.npz')==r['model_arrays_sha256'];d=f.data();own=np.array([list(d['genes']).index(t) for t in d['targets']])
 with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as z:
  assert z['effects'].shape==(3,300,18533);assert np.isfinite(z['effects']).all();assert (z['effects'][:,np.arange(300),own]==0).all();correction=float(abs(z['effects']-z['parent_effects']).max());assert correction<=.125002
 for c,tr in r['guard_traces'].items():assert tr['v66']['all_target_margin_constraints_passed'] and tr['v90']['all_target_margin_constraints_passed']
 write(A/'model_audit.json',dict(passed=True,gaussian_actual_source_parity_cases=len(ids),gaussian_absolute_tolerance=1e-10,priors=priors,all_effects_finite=True,all_own_transfer_effects_zero=True,own_factor_applied_once_in_emitter=.1,maximum_interval_correction=correction,identity_margin_guards_passed=True,source_cache_sha256=sha(A/'source_cache.npz'),model_sha256=sha(A/'reference_model.json'),final_accuracy_claim=False,auditor_sha256=sha(Path(__file__))));print('Final response numerical/provenance audit passed',flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
