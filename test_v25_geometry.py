"""Test geometry invariants and falsify covariance-only universal response identification."""
import unittest,numpy as np
from v25_methods import operator,covariance,force_confidence
class Geometry(unittest.TestCase):
 def test_identity_and_bures_map(self):
  s=np.diag([1.,2.,3.]);q=np.diag([1.2,2.2,3.2]);np.testing.assert_allclose(operator(s,s,'force'),np.eye(3),atol=1e-12);a=operator(q,s,'bures');np.testing.assert_allclose(a@s@a.T,q,atol=1e-10)
 def test_bounds_and_split_abstention(self):
  c=covariance(np.random.default_rng(25).normal(size=(100,5)));a=operator(c*100,c,'force');v=np.linalg.svd(a,compute_uv=False);self.assertLessEqual(v.max(),1.500001);np.testing.assert_array_equal(force_confidence(np.ones(3),-np.ones(3)),np.zeros(3))
 def test_same_geometry_different_intervention(self):
  force=np.array([-1.,0.]);responses=[]
  for omega in [0.,2.]:
   drift=np.array([[-1.,-omega],[omega,-1.]]);sigma=np.eye(2);np.testing.assert_array_equal(drift@sigma+sigma@drift.T+2*np.eye(2),np.zeros((2,2)));responses.append(-np.linalg.solve(drift,force))
  self.assertGreater(np.linalg.norm(responses[0]-responses[1]),.5)
if __name__=='__main__':unittest.main()
