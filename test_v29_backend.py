"""Cross-backend exact integer parity over varied feasible margin problems."""
import unittest
import numpy as np
from v29_pairwise_repair import repair_pairwise
from v29_pairwise_reference import repair_pairwise as reference
from v26_methods import integerize
class BackendParity(unittest.TestCase):
 def test_exact_parity(self):
  for seed in range(32):
   rng=np.random.default_rng(seed);n=12+seed%9;g=17+seed%11;ref=rng.poisson(rng.uniform(1,30,(n,1)),(n,g)).astype(np.int32);q=rng.lognormal(size=ref.shape);q/=q.sum(1,keepdims=True);x,_=integerize(q,ref.sum(1),ref.sum(0),seed);a,d=repair_pairwise(x,ref,max_fraction=.1);bb,e=reference(x,ref,max_fraction=.1);np.testing.assert_array_equal(a,bb);self.assertEqual(d['CPM_repair_molecules'],e['CPM_repair_molecules']);self.assertEqual(d['CPM_repair_steps'],e['CPM_repair_steps'])
if __name__=='__main__':unittest.main()
