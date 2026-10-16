"""Check that post-mixture calibration applies own-gene knockdown exactly once."""
import unittest
from unittest.mock import patch
import numpy as np
from scipy import sparse
import v37_recalibrated as m

class RecalibrationTests(unittest.TestCase):
 def check_emitter(self,m):
  template=sparse.csr_matrix([[8,2]]*4);libs=np.full(4,10);p=np.array([.8,.2]);zero=np.zeros(2)
  class Atlas:
   def dual_moment_counts(self,initial,wp,wb,depths,seed,iterations):
    np.testing.assert_allclose(initial.sum(1),1)
    np.testing.assert_allclose(wp,[2/7,5/7]);np.testing.assert_allclose(wb,wp)
    np.testing.assert_array_equal(depths,libs)
    assert seed==17 and iterations==200
    return np.array([[3,7]]*4)
  with patch.object(m,'original',return_value=(template,{'mean_probability_L1':0.,'bulk_probability_L1':0.})):
   out,diag=m.emit(template,zero,zero,zero,0,17,(libs,p,p),Atlas(),m.NAME)
  np.testing.assert_array_equal(out.sum(1).A1,libs)
  self.assertEqual(out.dtype.kind,'i');self.assertTrue((out.data>=0).all())

 def test_population_recalibration(self):
  self.check_emitter(m)
 def test_noise_recalibration(self):
  import v37_noise_recalibrated
  self.check_emitter(v37_noise_recalibrated)

if __name__=='__main__':unittest.main()
