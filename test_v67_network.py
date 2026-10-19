"""Scientific invariants for V67's new intervention-network operator."""
import unittest
import numpy as np
from threadpoolctl import threadpool_limits
from v67_intervention_network import graph, susceptibility, response_map, availability, correction


class NetworkTests(unittest.TestCase):
    def test_identical_context_is_exact_identity(self):
        w=np.array([[0,.2,0],[0,0,-.1],[.1,0,0.]])
        np.testing.assert_array_equal(response_map(w,np.ones(3)),np.eye(3))
        p=np.array([10.,20.,100.])
        np.testing.assert_array_equal(availability(p,p),np.ones(3))

    def test_zero_network_cannot_invent_response(self):
        np.testing.assert_array_equal(response_map(np.zeros((3,3)),np.array([.5,1.,1.5])),np.eye(3))

    def test_recover_known_feedback_network(self):
        w=np.array([[0,.1,0],[.12,0,-.1],[.05,0,0.]])
        s=susceptibility(w,np.ones(3))
        r=(s-np.eye(3)).T*np.log2(.1)
        got,_=graph(r,np.zeros_like(r),np.ones_like(r,bool))
        np.testing.assert_allclose(got,w,atol=1e-14)
        q=np.array([.7,1.2,1.1])
        np.testing.assert_allclose(s@response_map(got,q),susceptibility(w,q),atol=1e-14)

    def test_mediator_chain_propagates_context_change(self):
        w=np.array([[0,0,0],[.2,0,0],[0,.3,0]])
        hs=susceptibility(w,np.ones(3));hq=susceptibility(w,np.array([1.,.5,1.]))
        self.assertAlmostEqual(hq[2,0],hs[2,0]*.5)
        self.assertAlmostEqual(hq[1,0],hs[1,0])

    def test_stable_under_bounded_query_gains(self):
        rng=np.random.default_rng(67067);r=rng.normal(size=(30,30))
        w,trace=graph(r,np.zeros_like(r),np.ones_like(r,bool))
        self.assertLessEqual(trace['bounded_operator_norm'],.450000000001)
        h=susceptibility(w,rng.uniform(.5,1.5,30))
        self.assertTrue(np.isfinite(h).all())
        np.testing.assert_allclose(np.diag(h),1,atol=1e-14)

    def test_missing_pairs_abstain_and_true_zero_is_valid(self):
        a=availability(np.array([np.nan,0,30.,30.]),np.array([20.,20.,np.nan,0.]))
        np.testing.assert_array_equal(a,np.array([1.,.5,1.,1.]))

    def test_invalid_source_values_cannot_change_graph(self):
        r=np.array([[0,.2,.1],[.1,0,-.2],[0,.1,0.]])
        valid=np.ones_like(r,bool);valid[0,2]=False
        w,_=graph(r,np.zeros_like(r),valid)
        r[0,2]=1e6
        got,_=graph(r,np.zeros_like(r),valid)
        np.testing.assert_array_equal(w,got)


if __name__=='__main__':
    with threadpool_limits(limits=2):unittest.main()
