"""V78 exact row/column and zero-shape invariants on a varied-depth fixture."""
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
from v78_margin_variance import transform
from v69_count_bridge import parameters
rng=np.random.default_rng(7878);p=rng.dirichlet(np.full(128,.8));depth=rng.integers(2000,12000,64);x=sparse.csr_matrix(np.array([rng.multinomial(int(n),p) for n in depth],dtype=np.int32));q=parameters(x)
with threadpool_limits(limits=2):
 identical,diag=transform(x,q,np.zeros(128),'control10',7878)
 assert (identical!=x).nnz==0 and diag['zero_shape_identity']
 q=dict(q,alpha=np.clip(q['alpha']*1.5,.005,20))
 changed,diag=transform(x,q,np.zeros(128),'control25',7878)
 np.testing.assert_array_equal(changed.sum(0),x.sum(0));np.testing.assert_array_equal(changed.sum(1),x.sum(1));assert changed.dtype.kind in 'iu' and changed.data.min(initial=0)>=0 and diag['CPM_L1']<=.0005
 print('V78 exact count margins, zero-shape identity, CPM guard and raw integer checks passed.',diag,flush=True)

with threadpool_limits(limits=2):
 for arm in ('query50','query100'):
  out,dd=transform(x,q,np.zeros(128),arm,7878)
  np.testing.assert_array_equal(out.sum(0),x.sum(0));np.testing.assert_array_equal(out.sum(1),x.sum(1));assert dd['CPM_L1']<=.0005
 print('V78 full query-shape candidates also preserve exact margins and CPM safeguards.')
