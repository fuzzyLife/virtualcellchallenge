"""Focused checks for V5 ratio uncertainty, distinct moments, leakage and integer count generation."""
import unittest
import numpy as np
from scipy import sparse
from concurrent.futures import ThreadPoolExecutor
from threadpoolctl import threadpool_limits
from fit_v5_submission import moments
from build_v5_submission import desired_profiles,project_bulk,emit
from calibrate_v5_promoters import predictions,ramp
from evaluate_atlas_transfer import load_atlas_model


class V5Tests(unittest.TestCase):
    def test_pooled_uncertainty_accounts_for_depth_covariance(self):
        x=sparse.csr_matrix([[10,30],[20,60],[40,120]],dtype=float)
        probability,variance=moments(x)
        np.testing.assert_allclose(probability,[.25,.75])
        np.testing.assert_allclose(variance,0,atol=1e-18)

    def test_mean_and_pooled_profiles_remain_distinct_and_knockdown_is_relative(self):
        mean=np.array([.2,.3,.5]);bulk=np.array([.4,.1,.5]);zero=np.zeros(3)
        p,b=desired_profiles(mean,bulk,zero,zero,0)
        self.assertAlmostEqual(p[0]/p[1],.1*mean[0]/mean[1])
        self.assertAlmostEqual(b[0]/b[1],.1*bulk[0]/bulk[1])
        self.assertGreater(np.abs(p-b).sum(),.1)
        self.assertAlmostEqual(p.sum(),1);self.assertAlmostEqual(b.sum(),1)
        np.testing.assert_array_equal(mean,[.2,.3,.5])

    def test_promoter_strength_excludes_held_targets_and_preserves_stronger_response(self):
        left=np.array([-1.,-2.,100.,200.]);right=left-.2;scale=np.ones(4);train=np.array([True,True,False,False])
        pl,pr=predictions(left,right,scale,train,1.)
        left[~train]*=-1000;right[~train]*=1000
        a,b=predictions(left,right,scale,train,1.)
        np.testing.assert_array_equal(a,pl);np.testing.assert_array_equal(b,pr)
        a,b=predictions(left,right,scale,train,1.,np.repeat(-3.,4),np.zeros(4))
        np.testing.assert_array_equal(a,np.repeat(-3.,4));self.assertTrue((b<0).all())
        self.assertEqual(ramp(5000),0);self.assertEqual(ramp(100),1)

    def test_integer_emission_preserves_depth_and_is_parallel_deterministic(self):
        model,_=load_atlas_model();rng=np.random.default_rng(51);depth=np.linspace(10000,30000,40).astype(int)
        probability=rng.dirichlet(np.ones(12),40);x=sparse.csr_matrix(np.stack([rng.multinomial(int(n),p) for n,p in zip(depth,probability)]))
        p=np.mean(x.toarray()/depth[:,None],axis=0);b=np.asarray(x.sum(axis=0)).ravel()/depth.sum();rows=np.arange(40)
        with threadpool_limits(limits=2):
            baseline,_=emit(model,x,depth,rows,p,b,777)
            with ThreadPoolExecutor(max_workers=2) as pool:
                results=list(pool.map(lambda _:emit(model,x,depth,rows,p,b,777)[0],range(2)))
        for result in results:
            self.assertEqual((result-baseline).nnz,0)
            np.testing.assert_array_equal(np.asarray(result.sum(axis=1)).ravel(),depth)
        self.assertTrue(np.issubdtype(baseline.dtype,np.integer) and (baseline.data>0).all())

    def test_equal_depth_incompatible_profiles_are_rejected(self):
        p=np.array([.5,.5]);b=np.array([.9,.1]);depth=np.array([10000,10000])
        np.testing.assert_allclose(project_bulk(p,b,depth),p)
        model,_=load_atlas_model()
        with self.assertRaisesRegex(ValueError,'projection'):
            emit(model,sparse.csr_matrix([[5000,5000],[5000,5000]]),depth,np.arange(2),p,b,1)


if __name__=='__main__':unittest.main()
