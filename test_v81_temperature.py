"""Exact-count-margin, deterministic replay and zero-residual checks for V81."""
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v81_residual_temperature as m
from v57_pie import write,sha
rng=np.random.default_rng(81);p=rng.dirichlet(np.ones(128));x=sparse.csr_matrix(rng.poisson(rng.integers(1000,8000,40)[:,None]*p),dtype=np.int32)
with threadpool_limits(limits=2):
 for arm in m.ARMS:
  for opposite in [False,True]:
   y,d=m.transform(x,None,None,arm,81,opposite);z,_=m.transform(x,None,None,arm,81,opposite)
   assert d['CPM_L1']<=.0005 and y.dtype==np.int32 and (y.data>=0).all()
   np.testing.assert_array_equal(y.sum(0),x.sum(0));np.testing.assert_array_equal(y.sum(1),x.sum(1));assert (y!=z).nnz==0
 balanced=sparse.csr_matrix(np.arange(1,7)[:,None]*np.array([[0,2,4,6,8,10]],np.int32))
 for arm in m.ARMS:
  y,_=m.transform(balanced,None,None,arm,81);assert (y!=balanced).nnz==0
write(m.A/'additional_numerical_tests.json',dict(passed=True,deterministic_replay=True,balanced_zero_residual_identity=True,exact_integer_row_and_column_margins=True,CPM_bound=True,test_sha256=sha(__file__)))
print('V81 margins, zero-residual identity and deterministic replay passed',flush=True)
