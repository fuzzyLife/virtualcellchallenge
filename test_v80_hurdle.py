"""Hurdle-control fitting and exact-profile projection invariants."""
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
from v80_hurdle_controls import fit_hurdle,transform
from v69_count_bridge import parameters
rng=np.random.default_rng(8080);p=rng.dirichlet(np.full(96,.6));depth=rng.integers(2000,12000,64);x=sparse.csr_matrix(np.array([rng.multinomial(int(n),p) for n in depth],dtype=np.int32));q=parameters(x);q['genes']=np.asarray([str(i) for i in range(96)])
fit=fit_hurdle(x,q);assert np.isfinite(fit['hazard']).all() and (fit['alpha_positive']>=0).all()
with threadpool_limits(limits=2):
 for arm in ('hurdle_coupled','hurdle_independent','empirical'):
  out,diag=transform(x,fit,np.zeros(96),arm,8080,template=x)
  np.testing.assert_array_equal(out.sum(0),x.sum(0));np.testing.assert_array_equal(out.sum(1),x.sum(1));assert diag['CPM_L1']<=.0005
print('V80 empirical detection/positive-shape fitting and exact count/profile safeguards passed.',flush=True)
