"""Numerical tests for the independently implemented count bridge."""
import unittest
import numpy as np
from scipy import sparse
import v69_count_bridge as m

class TestBridge(unittest.TestCase):
    def test_null_exact(self):
        x=np.random.default_rng(1).poisson(3,size=(120,12));p=m.parameters(x)
        y,diag=m.bridge(x,np.zeros(12),None,9,p)
        np.testing.assert_array_equal(y.toarray(),x);self.assertTrue(diag['zero_response_identity'])
    def test_reproducible_integer_and_no_forced_margin(self):
        x=np.random.default_rng(2).negative_binomial(2,.4,size=(80,8));p=m.parameters(x);e=np.arange(8)/10-.3
        a,_=m.bridge(x,e,2,42,p);b,_=m.bridge(x,e,2,42,p)
        np.testing.assert_array_equal(a.toarray(),b.toarray());self.assertEqual(a.dtype.kind,'i')
        self.assertTrue((a.data>=0).all());self.assertTrue(np.any(a.sum(1).A1!=x.sum(1)))
    def test_moment_estimate(self):
        rng=np.random.default_rng(3);x=rng.negative_binomial(4,.5,size=(40000,6));p=m.parameters(x)
        # Conditional library modelling absorbs some between-gene variation;
        # raw estimate must be finite, bounded, and empirical proportions exact.
        np.testing.assert_allclose(p['p'],x.sum(0)/x.sum())
        self.assertTrue(np.isfinite(p['alpha']).all());self.assertTrue((p['alpha']>=m.ALPHA_MIN).all())
    def test_shape_coupling_marginals(self):
        # Nearly fixed external depth supplied by a large inert anchor gene.
        rng=np.random.default_rng(4);n=80000;r0=4.;mu=8.
        obs=rng.negative_binomial(r0,r0/(r0+mu),size=n)
        x=np.column_stack([obs,np.full(n,100000)])
        p=dict(p=np.array([mu/(100000+mu),100000/(100000+mu)]),alpha=np.array([1/r0,.005]))
        for ratio in [.5,2.]:
            delta=np.array([np.log(ratio),0.]);y,_=m.bridge(x,np.zeros(2),None,17,p,delta,1.)
            v=y[:,0].toarray().ravel();expected_variance=mu+mu*mu*(1/r0)*ratio
            self.assertAlmostEqual(v.mean()/mu,1.,delta=.025)
            self.assertAlmostEqual(v.var()/expected_variance,1.,delta=.06)
    def test_zero_measured_gene_abstains(self):
        x=np.column_stack([np.full(100,5),np.zeros(100,dtype=int)])
        p=m.parameters(x);y,_=m.bridge(x,np.array([0.,2.]),None,5,p,np.ones(2),.5)
        self.assertEqual(y[:,1].nnz,0)
    def test_own_thinning_effect(self):
        x=np.full((1000,10),10);p=m.parameters(x);y,_=m.bridge(x,np.zeros(10),0,6,p)
        self.assertLess(y[:,0].mean(),2.);self.assertGreater(y[:,1].mean(),9.)

if __name__=='__main__':unittest.main()
