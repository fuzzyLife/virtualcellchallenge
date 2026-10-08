"""Scientific invariants for source exclusion and population count transport."""
import unittest
import numpy as np
from scipy import sparse
from v10_methods import emit,response_bank
from v6_methods import moments
from audit_transfer_validation import pds

class V10Tests(unittest.TestCase):
 def test_integer_transport_depth_and_two_moments(self):
  rng=np.random.default_rng(41)
  d=np.linspace(.4,2,400)
  x=sparse.csr_matrix(rng.poisson(d[:,None]*[3,10,25,60,8,1]))
  eff=np.array([.2,-.4,.1,0.,.3,0.])
  for mode in ['compact','selective']:
   y,diag=emit(x,np.arange(400),eff,5,23,mode)
   np.testing.assert_array_equal(y.sum(axis=1),x.sum(axis=1))
   self.assertTrue(np.issubdtype(y.dtype,np.integer))
   self.assertTrue((y.data>0).all())
   self.assertLess(diag['mean_probability_L1'],.03)
   self.assertLess(diag['bulk_probability_L1'],.003)
   self.assertEqual((y-emit(x,np.arange(400),eff,5,23,mode)[0]).nnz,0)
 def test_zero_response_legacy_identity(self):
  x=sparse.csr_matrix([[1,2,3],[4,2,1]])
  y,_=emit(x,np.arange(2),np.zeros(3),None,12,'legacy')
  self.assertEqual((y-x).nnz,0)
 def test_population_not_control_vector_copy(self):
  rng=np.random.default_rng(8)
  delta=rng.integers(-2,3,size=400)
  x=sparse.csr_matrix(np.c_[3+delta,7-delta,np.full(400,15),np.full(400,50)])
  a,_=emit(x,np.arange(400),np.zeros(4),None,42,'compact')
  b,_=emit(x,np.arange(399,-1,-1),np.zeros(4),None,42,'compact')
  # Individual expression changes but equal depths and population moments do
  # not: the compact generator must not secretly copy original cell vectors.
  self.assertEqual((a-b).nnz,0)
  self.assertFalse((a-x).nnz==0)
 def test_held_family_responses_cannot_change_pool(self):
  from unittest.mock import patch
  families=np.array(['K562','CD4','H1','HCT116','HEK293T'])
  data={'targets':np.array(['a','b']), 'genes':np.array(['a','b','c','d']),
        'row_families':np.repeat(families,2),'row_targets':np.tile(np.arange(2),5),
        'truth':np.ones((10,4))*.2,'errors':np.ones((10,4))*.1,
        'valid':np.ones((10,4),bool),'knockdown':np.ones(10)*.7}
  with patch('v10_methods.old_response',return_value=np.zeros((2,4))):
   first,_=response_bank(data,excluded='HCT116')
   data['truth'][data['row_families']=='HCT116']=10000
   second,detail=response_bank(data,excluded='HCT116')
  self.assertNotIn('HCT116',detail['training_families'])
  for mode in ['pooled','dose','invariant']:np.testing.assert_array_equal(first[mode],second[mode])
 def test_constant_profile_has_chance_discrimination(self):
  t=np.array([[1.,2,3],[2.,1,3],[3.,2,1]])
  self.assertAlmostEqual(pds(np.ones_like(t),t,np.ones(3)*.5),.5)
if __name__=='__main__':unittest.main(verbosity=2)
