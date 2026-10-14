"""Scientific tests of V23 integer mean preservation and abstention."""
import unittest,numpy as np
from scipy import sparse
from v23_mean_repair import preserve_reference_means
from v23_methods import emit
from v6_methods import moments
from evaluate_atlas_transfer import load_atlas_model
class Repair(unittest.TestCase):
 def test_integer_conservation_and_precision(self):
  rng=np.random.default_rng(23);ref=rng.poisson(10,size=(40,80)).astype(np.int32);x=ref.copy()
  for i in range(500):
   row=int(rng.integers(40));a,b=rng.choice(80,2,replace=False)
   if x[row,a]>0:x[row,a]-=1;x[row,b]+=1
  a,d=preserve_reference_means(x,ref);b,e=preserve_reference_means(x,ref);np.testing.assert_array_equal(a,b);np.testing.assert_array_equal(a.sum(0),ref.sum(0));np.testing.assert_array_equal(a.sum(1),ref.sum(1));self.assertLess(d['mean_CPM_probability_L1_after'],d['mean_CPM_probability_L1_before']/100)
  same,_=preserve_reference_means(ref,ref);np.testing.assert_array_equal(same,ref)
  with self.assertRaises(AssertionError):preserve_reference_means(x.astype(float)+.5,ref)
 def test_zero_response_exact_baseline(self):
  x=sparse.csr_matrix(np.random.default_rng(24).poisson(30,size=(40,16)));libs,p,b,_=moments(x);atlas,_=load_atlas_model();base,_=emit(x,np.zeros(16),None,24,'v16',(libs,p,b),atlas);new,d=emit(x,np.zeros(16),None,24,'mean_preserved',(libs,p,b),atlas);self.assertEqual((base!=new).nnz,0)
if __name__=='__main__':unittest.main()
