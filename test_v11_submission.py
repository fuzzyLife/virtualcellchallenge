"""Scientific V11 tests: V10 equivalence, exclusion, consensus and depth integrity."""
import unittest
from unittest.mock import patch
import numpy as np
from scipy import sparse
from v11_methods import emit,response_bank,consensus_gate,SETTINGS
from v10_methods import emit as emit10,response_bank as bank10
from v6_methods import moments

class V11Tests(unittest.TestCase):
 def test_v10_emission_exact_and_nondefault_constraints(self):
  rng=np.random.default_rng(41);x=sparse.csr_matrix(rng.poisson(np.linspace(.4,2,400)[:,None]*[3,10,25,60,8,1]))
  effect=np.array([.2,-.4,.1,0.,.3,0.])
  baseline,_=emit10(x,np.arange(400),effect,5,23,'compact')
  current,_=emit(x,np.arange(400),effect,5,23,SETTINGS['v10'])
  self.assertEqual((baseline-current).nnz,0)
  for name in ['noise35','noise65','bulk085']:
   count,diag=emit(x,np.arange(400),effect,5,23,SETTINGS[name])
   np.testing.assert_array_equal(count.sum(axis=1),x.sum(axis=1));self.assertTrue(np.issubdtype(count.dtype,np.integer))
   self.assertTrue((count.data>0).all());self.assertLess(diag['bulk_probability_L1'],.003)
   self.assertEqual((count-emit(x,np.arange(400),effect,5,23,SETTINGS[name])[0]).nnz,0)
 def test_wrong_template_depth_adapter_rejected(self):
  x=sparse.csr_matrix(np.ones((400,4),dtype=int))
  with self.assertRaises(ValueError):emit(x,np.arange(400),np.zeros(4),None,3,SETTINGS['v10'],prepared=(np.ones(800)*4,np.ones(4)/4,np.ones(4)/4))
 def test_consensus_needs_independent_direction_support(self):
  y=np.array([[[1.,1.,1.]],[[1.,-1.,1.]]]);valid=np.ones_like(y,bool);valid[1,0,2]=False
  gate,_,_=consensus_gate(y,np.ones_like(y)*.1,valid,np.ones((2,1)),np.ones(2),np.ones((1,3)),np.ones((1,3)))
  self.assertGreater(gate[0,0],0);self.assertEqual(gate[0,1],0);self.assertEqual(gate[0,2],0)
 def test_baseline_response_and_held_family_invariance(self):
  families=np.array(['K562','CD4','H1','HCT116','HEK293T']);rng=np.random.default_rng(31)
  data={'targets':np.array(['a','b']),'genes':np.array(['a','b','c','d']),
   'row_families':np.repeat(families,2),'row_targets':np.tile(np.arange(2),5),'truth':rng.normal(0,.2,(10,4)),
   'errors':np.ones((10,4))*.1,'valid':np.ones((10,4),bool),'knockdown':np.ones(10)*.7}
  with patch('v11_methods.old_response',return_value=np.zeros((2,4))),patch('v10_methods.old_response',return_value=np.zeros((2,4))):
   first,_=response_bank(data,excluded='HCT116');old,_=bank10(data,excluded='HCT116')
   np.testing.assert_array_equal(first['dose'],old['dose'])
   data['truth'][data['row_families']=='HCT116']=10000
   second,record=response_bank(data,excluded='HCT116');self.assertNotIn('HCT116',record['training_families'])
   for key in ['dose','consensus20','consensus40']:np.testing.assert_array_equal(first[key],second[key])
if __name__=='__main__':unittest.main(verbosity=2)
