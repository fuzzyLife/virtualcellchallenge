"""Scientific invariants for the corrected V8 evaluation and count decoder."""
import unittest
import numpy as np,pandas as pd
from scipy import sparse
from audit_transfer_validation import self_test,pds
from prepare_v8_validation import split_controls,global_target_split
from v8_methods import emit,gated
from v6_methods import moments,direction_metrics
class V8Tests(unittest.TestCase):
 def test_official_profile_limits(self):self.assertEqual(self_test()['constant_prediction_pds'],.5)
 def test_shared_control_required(self):
  with self.assertRaises(ValueError):pds(np.ones((4,8)),np.ones((4,8)),np.ones((4,8)))
 def test_control_split_does_not_read_target_mixture(self):
  m=pd.DataFrame({'sample':['a']*16+['b']*16+['a']*8,'gene_target':['Non-Targeting']*32+['target']*8})
  left,right=split_controls(m,11);m.loc[32:,'sample']='b'
  a,b=split_controls(m,11);np.testing.assert_array_equal(left,a);np.testing.assert_array_equal(right,b)
  self.assertFalse(set(left)&set(right));self.assertEqual(len(left),16)
 def test_target_assignment_stable_across_contexts(self):
  self.assertEqual(global_target_split('VIM'),global_target_split('VIM'))
 def test_gating_excludes_own(self):
  result=gated(np.array([1.,2.,3.]),np.array([10.,10.,0.]),0)
  np.testing.assert_array_equal(result,[0.,2.,0.])
 def test_neutral_legacy_preserves_counts(self):
  x=sparse.csr_matrix([[0,3,7],[1,2,9]],dtype=np.int32)
  result,diag=emit(x,np.arange(2),np.zeros(3),None,19,'legacy')
  self.assertEqual((result!=x).nnz,0)
 def test_dual_preserves_depths_and_profile(self):
  rng=np.random.default_rng(7);x=sparse.csr_matrix(rng.poisson([4,10,30,8],size=(64,4)))
  for name in ['dual','dual_half']:
   result,diag=emit(x,np.arange(64),np.log2([1.1,.9,1.,1.]),None,23,name)
   np.testing.assert_array_equal(np.asarray(result.sum(axis=1)).ravel(),np.asarray(x.sum(axis=1)).ravel())
   self.assertTrue(np.issubdtype(result.dtype,np.integer));self.assertTrue((result.data>=0).all())
   self.assertLess(diag['projected_bulk_probability_L1'],.002)
 def test_safety_fidelity_uses_common_reference_budget(self):
  from select_v8_safe_generator import supported_budget_fidelity
  frame=pd.DataFrame({'candidate':['a','a','b','b'],'reference_calls':[100,1,100,1],
                      'predicted_calls':[100,1,100,1],'fidelity':[.4,1.,.4,0.]})
  self.assertEqual(supported_budget_fidelity(frame),{'a':.4,'b':.4})
 def test_fidelity_accounts_for_response_budget(self):
  genes=np.array(['own','a','b','c']);truth={'eligible':np.ones(4,bool),'q':np.array([0.,.01,.01,.9]),'lfc':np.array([-1.,1.,-1.,1.])}
  pred={'p':np.ones(4)*.001,'q':np.array([0.,.01,.9,.9]),'lfc':np.array([-1.,1.,1.,1.])}
  r=direction_metrics(pred,truth,genes,0);self.assertEqual(r['fidelity'],.5);self.assertEqual(r['predicted_calls'],1);self.assertEqual(r['reference_calls'],2)
if __name__=='__main__':unittest.main(verbosity=2)
