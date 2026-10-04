"""Tests for source-aware blending, uncertainty shrinkage and named contexts."""
import unittest
import numpy as np
from fit_context_ensemble import target_folds, estimate, combine_sources, CONDITIONS
from build_submission import context_factors

class ContextEnsembleTests(unittest.TestCase):
    def test_target_folds_do_not_split_conditions(self):
        targets=np.array(['T1','T2','T1','T3','T2','T4'])
        folds=target_folds(targets,3)
        self.assertEqual(folds[0],folds[2]); self.assertEqual(folds[1],folds[4])
        np.testing.assert_array_equal(folds[::-1],target_folds(targets[::-1],3))

    def test_uncertainty_shrinkage_and_zero_anchor(self):
        x=np.array([[1.,-1.,2.]])
        se=np.array([[.2,.2,2.]])
        q=np.array([[.01,.01,.9]])
        c=np.array(['Rest']); prior={'Rest':np.zeros(3)}
        settings={'kind':'eb','tau2':.04,'alpha':.5,'prior':False}
        actual=estimate(x,se,q,c,settings,prior)
        np.testing.assert_allclose(actual,[ [.25,-.25, .04/4.04] ])
        np.testing.assert_array_equal(estimate(x,se,q,c,{'kind':'zero'},prior),np.zeros_like(x))

    def test_blending_renormalizes_measured_gene_support(self):
        cd4=np.array([[[1.,1.,1.,1.]],[[2.,2.,2.,2.]],[[3.,3.,3.,3.]]])
        ca=np.ones((3,1),bool)
        k=np.full((1,4),4.); h=np.full((1,4),6.)
        mixture={context:{'cell_profiles':{'CD4':.6,'K562':.2,'H1':.2,'HCT116':0.,'HEK293T':0.},
            'cd4_conditions':dict(zip(['CD4_'+c for c in CONDITIONS],[.2,.3,.5]))} for context in ['A','B','C']}
        effects,detail=combine_sources(cd4,ca,k,np.array([True]),h,np.array([True]),mixture,
            np.array([True,False,False,False]),np.array([False,True,False,False]),np.array([[False,False,True,False]]))
        np.testing.assert_allclose(effects,np.tile([2.15,4.,6.,0.],(3,1,1)),atol=1e-6)
        self.assertEqual(len(detail),3)
        # A missing source condition must not dilute the available estimates with zero.
        ca[0,0]=False
        effects,_=combine_sources(cd4,ca,k,np.array([False]),h,np.array([False]),mixture)
        np.testing.assert_allclose(effects,np.full((3,1,4),113/44),atol=1e-6)

    def test_context_specific_model_requires_named_context(self):
        effects=np.array([[[1.,0.]],[[0.,1.]],[[.5,.5]]])
        summary={'context_specific':True,'model':'test','context_mixture_weights':{'A':{},'B':{},'C':{}}}
        model=(effects,np.ones((3,1),bool),np.zeros((3,2)),np.arange(2),{},['A','B','C'],summary)
        with self.assertRaises(ValueError):
            context_factors(np.ones(2),model,np.array(['G','T']),np.array(['T']),.9)
        a,_=context_factors(np.ones(2),model,np.array(['G','T']),np.array(['T']),.9,context='A')
        b,_=context_factors(np.ones(2),model,np.array(['G','T']),np.array(['T']),.9,context='B')
        np.testing.assert_allclose(a,[[2.,.1]])
        np.testing.assert_allclose(b,[[1.,.1]])

if __name__=='__main__':
    unittest.main()
