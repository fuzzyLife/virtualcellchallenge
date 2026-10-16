"""Boundary and meaningful population-mixture tests for the ensemble."""
import unittest
import numpy as np
from scipy import sparse
from unittest.mock import patch
from v37_ensemble import combine,emit,ARMS

class EnsembleTests(unittest.TestCase):
 def test_consensus_and_convexity(self):
  b=np.array([[0.,0.,2.9,-2.9]],np.float32);r=b+np.array([[.1,-.1,.05,-.05]],np.float32);p=b+np.array([[.2,.2,.02,-.02]],np.float32);c=b.copy()
  z=combine(b,r,p,c)
  np.testing.assert_allclose(z['consensus']-b,[[.1,0,.02,-.02]],atol=2e-7)
  np.testing.assert_allclose(z['mean50'],(b+r)/2,atol=2e-7)
  for v in z.values():self.assertTrue(np.isfinite(v).all());self.assertLessEqual(abs(v).max(),3)
 def test_identical_experts_fallback(self):
  b=np.array([[0,.1,-.2]],np.float32)
  for v in combine(b,b,b,b).values():np.testing.assert_array_equal(v,b)
 def test_population_mixture_preserves_rows_and_libraries(self):
  t=sparse.csr_matrix([[2,3],[4,1],[1,4],[3,2]])
  aa=t.copy();bb=sparse.csr_matrix([[5,0],[0,5],[4,1],[1,4]])
  def parent(*args,**kwargs):return (aa if float(args[1][1])==0 else bb),{}
  with patch('v36_methods.emit',side_effect=parent):
   args=(t,np.array([0.,.1]),np.array([0.,0.]),np.array([0.,.2]),0,17,(np.full(4,5),np.array([.5,.5]),np.array([.5,.5])),None,'population50')
   out,d=emit(*args);again,_=emit(*args)
  np.testing.assert_array_equal(out.sum(1).A1,[5]*4)
  self.assertEqual((out!=aa).getnnz(1).astype(bool).sum(),2)
  self.assertEqual((out!=again).nnz,0)
  self.assertEqual(d['mixture_rows'],2)

if __name__=='__main__':unittest.main()
