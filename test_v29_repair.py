"""Pair-aware swap monotonicity, conservation, determinism and zero response."""
import unittest
import numpy as np
from v26_methods import integerize
from v29_pairwise_repair import repair_pairwise
class Pairwise(unittest.TestCase):
 def test_identity(self):
  x=np.random.default_rng(29).poisson(4,(40,51)).astype(np.int32);y,d=repair_pairwise(x,x);np.testing.assert_array_equal(y,x);self.assertEqual(d['CPM_repair_steps'],0)
 def test_conservation_monotone_budget_deterministic(self):
  rng=np.random.default_rng(30);ref=rng.poisson(np.linspace(2,30,40)[:,None],(40,51)).astype(np.int32);q=rng.uniform(size=ref.shape);q/=q.sum(1,keepdims=True);x,_=integerize(q,ref.sum(1),ref.sum(0),30);a,d=repair_pairwise(x,ref,max_fraction=.05);bb,e=repair_pairwise(x,ref,max_fraction=.05);np.testing.assert_array_equal(a,bb);np.testing.assert_array_equal(a.sum(1),ref.sum(1));np.testing.assert_array_equal(a.sum(0),ref.sum(0));self.assertLessEqual(d['mean_reference_L1'],d['initial_CPM_L1']+1e-10);self.assertLessEqual(d['CPM_repair_molecules'],d['CPM_repair_budget'])
 def test_equal_depth_feasible_identity_weighted_means(self):
  ref=np.array([[5,4,1],[2,3,5]],np.int32);x=np.array([[3,4,3],[4,3,3]],np.int32);y,d=repair_pairwise(x,ref);np.testing.assert_array_equal(y,x);self.assertLess(d['mean_reference_L1'],1e-12)
if __name__=='__main__':unittest.main()
