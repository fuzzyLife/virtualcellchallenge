"""Check activation preservation, partial moment fits and neutral identity."""
import unittest
import numpy as np
from scipy import sparse
from v9_methods import emit,GENERATORS
from v8_methods import emit as raw_emit
from v6_methods import moments

class V9Tests(unittest.TestCase):
 def test_neutral_legacy_is_original(self):
  x=sparse.csr_matrix([[0,3,7],[1,2,9]],dtype=np.int32)
  y,_=emit(x,np.arange(2),np.zeros(3),None,19)
  self.assertEqual((x-y).nnz,0)
 def test_partial_fits_preserve_simulated_libraries(self):
  rng=np.random.default_rng(8);x=sparse.csr_matrix(rng.poisson([1,20,100,50],size=(400,4)))
  effect=np.log2([4,.8,1.,1.]);initial,_=emit(x,np.arange(400),effect,1,23)
  for g in GENERATORS[1:]:
   y,d=emit(x,np.arange(400),effect,1,23,g)
   np.testing.assert_array_equal(y.sum(axis=1),initial.sum(axis=1))
   self.assertTrue(np.issubdtype(y.dtype,np.integer));self.assertTrue((y.data>=0).all())
   self.assertLess(d['projected_bulk_probability_L1'],.002)
   # An activated sparse gene must remain detected well beyond untreated prevalence.
   self.assertGreater(y[:,0].getnnz(),x[:,0].getnnz()+40)
 def test_full_neutral_matches_control_based_fit(self):
  rng=np.random.default_rng(4);x=sparse.csr_matrix(rng.poisson([3,8,20,100],size=(128,4)))
  a,_=emit(x,np.arange(128),np.zeros(4),None,17,'legacy_fit100')
  b,_=raw_emit(x,np.arange(128),np.zeros(4),None,17,'dual')
  self.assertEqual((a-b).nnz,0)
if __name__=='__main__':unittest.main(verbosity=2)
