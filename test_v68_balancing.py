"""Tests for constrained real-cell balancing and probability-preserving sampling."""
import unittest
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
from v68_population_balancing import strata,project,uniform,pivotal,sample,edit_own,operators,fit,permute_within


class BalancingTests(unittest.TestCase):
    def test_capped_stratum_projection(self):
        group,quota=strata(np.arange(50)+20,20)
        p=project(np.random.default_rng(1).normal(size=50),group,quota,20)
        self.assertGreaterEqual(p.min(),0);self.assertLessEqual(p.max(),.05+1e-12)
        for b,k in enumerate(quota):self.assertAlmostEqual(p[group==b].sum(),k/20,12)
        self.assertGreaterEqual(1/(p@p),20-1e-9)

    def test_pivotal_fixed_size_and_unbiased_inclusion(self):
        pi=np.array([.1,.4,.6,.9]);mean=np.zeros(4)
        for seed in range(6000):
            idx=pivotal(pi,np.random.default_rng(seed));self.assertEqual(len(idx),2);mean[idx]+=1
        np.testing.assert_allclose(mean/6000,pi,atol=.018,rtol=0)

    def test_sampling_exact_quota_and_unique(self):
        group,quota=strata(np.arange(80)+20,40)
        p=project(np.random.default_rng(2).normal(size=80),group,quota,40)
        ids=sample(p,group,quota,40,4)
        self.assertEqual(len(np.unique(ids)),40)
        np.testing.assert_array_equal(np.bincount(group[ids],minlength=len(quota)),quota)

    def test_own_edit_preserves_every_other_count(self):
        x=sparse.csr_matrix(np.arange(120).reshape(12,10)+1)
        y=edit_own(x,3,7);diff=(x-y).tocoo()
        self.assertTrue((diff.col==3).all());self.assertTrue((y.data>=0).all())
        np.testing.assert_array_equal(edit_own(x,None,4).toarray(),x.toarray())

    def test_operator_adjoint_and_objective_gradient(self):
        x=sparse.csr_matrix(np.random.default_rng(3).poisson(3,(30,12)))
        target=np.ones(12)/12;axis=np.ones(12,bool);f,t=operators(x,target,axis)
        rng=np.random.default_rng(8);a=rng.normal(size=30);b=rng.normal(size=12)
        self.assertAlmostEqual(float(f(a)@b),float(a@t(b)),10)
        w=np.ones(30)/30;d=rng.normal(size=30);eps=1e-7
        numerical=(.5*(f(w+eps*d)@f(w+eps*d))-.5*(f(w-eps*d)@f(w-eps*d)))/(2*eps)
        self.assertAlmostEqual(numerical,float(t(f(w))@d),7)

    def test_optimizer_improves_known_population_mixture(self):
        # Two control expression states, both in each depth stratum.
        x=np.tile(np.array([[8,2,5,5],[2,8,5,5]]),(40,1));x=sparse.csr_matrix(x)
        group=np.zeros(80,int);quota=np.array([20]);target=np.array([.32,.18,.25,.25])
        p,r=fit(x,target,np.ones(4,bool),20,group,quota)
        self.assertLess(r['final_objective'],r['initial_objective']*.1)
        self.assertTrue(r['converged']);self.assertGreater(p[::2].sum(),.65)

    def test_zero_response_matches_uniform_weight_solution(self):
        x=sparse.csr_matrix(np.tile([2,4,6,8],(40,1)));target=np.array([.1,.2,.3,.4])
        group,quota=strata(x.sum(1).A1,20);p,r=fit(x,target,np.ones(4,bool),20,group,quota)
        np.testing.assert_allclose(p,uniform(group,quota,20),atol=1e-14)
        self.assertLess(r['duality_gap'],1e-12)

    def test_shuffle_preserves_stratum_weight_distribution(self):
        group,quota=strata(np.arange(80)+20,40);p=project(np.random.default_rng(4).normal(size=80),group,quota,40)
        q=permute_within(p,group)
        for b in range(len(quota)):np.testing.assert_array_equal(np.sort(p[group==b]),np.sort(q[group==b]))


if __name__=='__main__':
    with threadpool_limits(limits=2):unittest.main()
