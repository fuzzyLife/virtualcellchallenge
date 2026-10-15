"""Exercise the frozen numerical bridge with an explicitly blocked two-gene swap."""
import unittest
from unittest.mock import patch
import numpy as np
from v29_methods import project
from v29_pairwise_repair import repair_pairwise
class Bridge(unittest.TestCase):
 def test_budgeted_carrier_resolves_blocked_pair(self):
  x=np.array([[0,0,100],[4,0,196],[0,4,396],[0,0,800]],np.int32)
  ref=np.array([[4,0,96],[0,0,200],[0,0,400],[0,4,796]],np.int32)
  _,d=repair_pairwise(x,ref);self.assertGreater(d['mean_reference_L1'],.0005);self.assertLess(d['CPM_repair_molecules'],d['CPM_repair_budget'])
  # Isolate the repair stage while preserving its real integer input and all
  # actual swap code. Probability fitting/initial rounding are tested elsewhere.
  with patch('v29_methods.fit_probabilities',return_value=x/x.sum(1)[:,None]),patch('v29_methods.integerize',return_value=(x.copy(),0)):
   y,diag=project(x/x.sum(1)[:,None],ref,29)
  self.assertEqual(diag['repair_path'],'joint_carrier_joint');self.assertLess(diag['mean_reference_L1'],.0005);np.testing.assert_array_equal(y.sum(1).A1,ref.sum(1));np.testing.assert_array_equal(y.sum(0).A1,ref.sum(0));self.assertLessEqual(diag['CPM_repair_molecules'],diag['CPM_repair_budget'])
if __name__=='__main__':unittest.main()
