"""Check guide-pooling invariants independently of real treatment outcomes."""
import unittest
import numpy as np
from v31_methods import teacher
class GuideTests(unittest.TestCase):
 def blocks(self):
  return [{'effect':np.array([1.,2.,3.]),'variance':np.array([.1,.2,.3]),'valid':np.array([True,True,False]),'n':80,'knockdown':.9},{'effect':np.array([2.,-1.,1000.]),'variance':np.array([.2,.4,.1]),'valid':np.array([True,False,False]),'n':120,'knockdown':.9}]
 def test_guide_order_invariance(self):
  a=self.blocks();x,c=teacher(a,'guide');y,d=teacher(a[::-1],'guide');np.testing.assert_allclose(x,y,rtol=1e-7);np.testing.assert_allclose(c,d,rtol=1e-7)
 def test_missing_coordinates_abstain(self):
  a=self.blocks();x,c=teacher(a,'guide');self.assertEqual(x[2],0);self.assertEqual(c[2],0);self.assertAlmostEqual(float(x[1]),2.)
 def test_unknown_knockdown_keeps_single_guide_mean(self):
  a=self.blocks()[:1];a[0]['knockdown']=None;x,_=teacher(a,'guide');np.testing.assert_array_equal(x[:2],a[0]['effect'][:2])
if __name__=='__main__':unittest.main()
