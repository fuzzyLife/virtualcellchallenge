"""Scientific invariants for the V22 conditional population proposal."""
import unittest
import numpy as np
from scipy import sparse
from v22_methods import conditional_delta,emit
from v16_methods import emit as old_emit
from v6_methods import moments
from evaluate_atlas_transfer import load_atlas_model

class StateResponse(unittest.TestCase):
 def model(self,g=16):
  return {'support':np.array([200]),'output_axis':np.array([2,3,4]),'strength_residual':np.array([np.linspace(-.2,.2,128)]),'strength_beta':np.array([[.1,-.1,0,0]]),'interaction':np.array([[[.2,0,.1],[0,.2,0],[.1,0,0],[0,0,.1]]])}
 def test_conditional_pairing(self):
  z=np.random.default_rng(7).normal(size=(400,4));effect=np.linspace(-.2,.2,16);m=self.model();q=np.full(16,20.)
  a=conditional_delta(z,effect,m,0,'conditional',11,q);b=conditional_delta(z,effect,m,0,'conditional',11,q)
  np.testing.assert_array_equal(a,b);self.assertLessEqual(abs(a).max(),.125)
  shuffled=conditional_delta(z,effect,m,0,'shuffled_state',11,q);self.assertGreater(np.linalg.norm(a-shuffled),.01)
  zero=dict(m,strength_residual=np.zeros((1,128)),strength_beta=np.zeros((1,4)),interaction=np.zeros((1,4,3)))
  np.testing.assert_array_equal(conditional_delta(z,effect,zero,0,'conditional',11,q),np.zeros((400,16)))
 def test_generator_null_and_libraries(self):
  rng=np.random.default_rng(8);x=sparse.csr_matrix(rng.poisson(50,size=(40,16)));libs,p,b,_=moments(x);atlas,_=load_atlas_model();model=self.model();z=rng.normal(size=(40,4))
  a,_=emit(x,np.zeros(16),None,9,'conditional',(libs,p,b),atlas,model,0,z);expected,_=old_emit(x,np.arange(40),np.zeros(16),None,9,'residual10',(libs,p,b),atlas)
  self.assertEqual((a!=expected).nnz,0)
  count,d=emit(x,np.linspace(-.2,.2,16),1,9,'conditional',(libs,p,b),atlas,model,0,z)
  np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),libs);self.assertTrue(np.issubdtype(count.dtype,np.integer));self.assertTrue((count.data>=0).all())
  abstain=dict(model,support=np.array([0]));a,_=emit(x,np.linspace(-.2,.2,16),1,9,'conditional',(libs,p,b),atlas,abstain,0,z);b,_=old_emit(x,np.arange(40),np.linspace(-.2,.2,16),1,9,'residual10',(libs,p,b),atlas);self.assertEqual((a!=b).nnz,0)
if __name__=='__main__':unittest.main()
