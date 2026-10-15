"""V27 numerical invariants, independent of biological metric outcomes."""
import unittest
import numpy as np
from scipy import sparse
from v27_methods import project,anchored_initial,relative_delta,control_calibration
class Invariants(unittest.TestCase):
 def test_identity_projection(self):
  rng=np.random.default_rng(27); x=rng.poisson(3,(24,31)).astype(np.int32); x[:,0]=0
  libs=x.sum(1); p=(x/libs[:,None]).mean(0)
  q=anchored_initial(sparse.csr_matrix(x),np.zeros(31),None,p)
  y,d=project(q,sparse.csr_matrix(x),27)
  np.testing.assert_array_equal(y.toarray(),x);self.assertLess(d['mean_reference_L1'],1e-12)
 def test_contract_and_determinism(self):
  rng=np.random.default_rng(28); ref=rng.poisson(4,(32,37)).astype(np.int32); q=rng.lognormal(size=ref.shape)
  x,d=project(q,ref,28);y,_=project(q,ref,28)
  np.testing.assert_array_equal(x.toarray(),y.toarray());np.testing.assert_array_equal(x.sum(1).A1,ref.sum(1));np.testing.assert_array_equal(x.sum(0).A1,ref.sum(0));self.assertTrue((x.data>=0).all())
 def test_quantile_calibration_deterministic_no_test_input(self):
  rng=np.random.default_rng(29); x=sparse.csr_matrix(rng.poisson(2,(20,31)));fit=sparse.csr_matrix(rng.poisson(3,(25,31)));p=np.ones(31)/31
  a=control_calibration(x,fit,p,29);bb=control_calibration(x,fit,p,29)
  np.testing.assert_array_equal(a,bb);np.testing.assert_allclose(a.sum(1),1);self.assertTrue((a>=0).all())
 def test_zero_condition_relative_identity(self):
  rng=np.random.default_rng(30); w={}
  for n,i,o in [('d1',40,64),('d2',64,128),('out',128,31)]:w[n+'.weight']=rng.normal(0,.01,(o,i)).astype(np.float32);w[n+'.bias']=np.zeros(o,np.float32)
  geo={'positions':np.arange(31),'scale':np.ones(31),'basis':rng.normal(size=(31,64))}
  model={'weights':w,'output_mask':np.ones(31,bool),'geometry':geo}
  delta=relative_delta(np.zeros((20,8),np.float32),np.zeros(31),np.ones(31)/31,model,30)
  np.testing.assert_array_equal(delta,np.zeros_like(delta))
if __name__=='__main__':unittest.main()
