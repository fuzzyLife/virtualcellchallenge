"""V28 deterministic global-only integer and mean contracts."""
import unittest
import numpy as np
from scipy import sparse
from v28_methods import project
class GlobalTransport(unittest.TestCase):
 def test_identity(self):
  r=np.random.default_rng(28).poisson(4,(40,51)).astype(np.int32);r[:,0]=0;libs=r.sum(1);q=r/libs[:,None];x,d=project(q,r,28);np.testing.assert_array_equal(x.toarray(),r);self.assertLess(d['mean_reference_L1'],1e-12)
 def test_exact_margins_determinism_and_no_depth_constraints(self):
  rng=np.random.default_rng(29);r=rng.poisson(np.linspace(2,30,40)[:,None],(40,51)).astype(np.int32);q=rng.uniform(size=r.shape);x,d=project(q,r,29);y,_=project(q,r,29);np.testing.assert_array_equal(x.toarray(),y.toarray());np.testing.assert_array_equal(x.sum(1).A1,r.sum(1));np.testing.assert_array_equal(x.sum(0).A1,r.sum(0));self.assertTrue((x.data>=0).all());self.assertFalse(d['depth_stratum_repair'])
if __name__=='__main__':unittest.main()
