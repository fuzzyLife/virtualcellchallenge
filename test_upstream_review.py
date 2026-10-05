"""Meaningful leakage and feasible-moment checks for the upstream review experiments."""
import unittest
import numpy as np
from evaluate_atlas_transfer import source_center,load_atlas_model


class UpstreamReviewTests(unittest.TestCase):
    def test_all_held_target_conditions_excluded_from_center(self):
        values=np.array([[1.,2.],[3.,4.],[100.,200.],[300.,400.]])
        conditions=np.array(['Rest','Stim8hr','Rest','Stim8hr'])
        valid=np.ones_like(values,dtype=bool);train=np.array([True,True,False,False])
        before=source_center(values,valid,conditions,train)
        values[~train]*=-1000
        np.testing.assert_array_equal(before,source_center(values,valid,conditions,train))
        np.testing.assert_array_equal(before,[[1,2],[3,4],[1,2],[3,4]])

    def test_own_or_unmeasured_gene_does_not_enter_center(self):
        values=np.array([[999.,2.],[3.,np.nan]])
        valid=np.array([[False,True],[True,False]])
        actual=source_center(values,valid,np.repeat('Rest',2),np.ones(2,bool))
        np.testing.assert_array_equal(actual,[[3,2],[3,2]])

    def test_integer_dual_moments_preserve_depth_and_pooled_profile(self):
        model,_=load_atlas_model();rng=np.random.default_rng(902)
        template=rng.gamma(2,1,size=(40,12));template/=template.sum(axis=1,keepdims=True)
        depths=np.arange(1000,1800,20);p=template.mean(axis=0)
        bulk=(depths[:,None]*template).sum(axis=0)/depths.sum()
        actual=model.dual_moment_counts(template,p,bulk,depths=depths,seed=10)
        repeat=model.dual_moment_counts(template,p,bulk,depths=depths,seed=10)
        np.testing.assert_array_equal(actual,repeat)
        np.testing.assert_array_equal(actual.sum(axis=1),depths)
        self.assertTrue(np.issubdtype(actual.dtype,np.integer) and (actual>=0).all())
        np.testing.assert_allclose(actual.sum(axis=0)/actual.sum(),bulk,rtol=0,atol=12/actual.sum())
        self.assertLess(np.abs((actual/depths[:,None]).mean(axis=0)-p).sum(),.003)

    def test_infeasible_equal_depth_moments_fail(self):
        model,_=load_atlas_model()
        with self.assertRaises(ValueError):
            model.dual_moment_counts([[10,0],[0,10]],[.5,.5],[.9,.1],depths=[100,100],seed=3)


if __name__=='__main__':unittest.main()
