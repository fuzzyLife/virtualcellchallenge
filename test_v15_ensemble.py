"""Scientific invariants for V15 ensemble and marginal transport."""
import unittest
import numpy as np
from scipy import sparse
from v15_ensemble_methods import cube,pool,source_reliability,GRID
from v15_population_methods import transport_initial,emit
from v6_methods import moments
from test_v14_context import ContextTests

class EnsembleTests(unittest.TestCase):
    def test_held_family_labels_are_excluded(self):
        d=ContextTests().fixture();a=cube(d,['HCT116']);p={k:v.copy() for k,v in d.items()};rows=p['row_families']=='HCT116'
        p['truth'][rows]=1e5;p['errors'][rows]=1e-8;p['knockdown'][rows]=.001;p['valid'][rows]=False
        b=cube(p,['HCT116'])
        for x,y in zip(a,b):np.testing.assert_array_equal(x,y)
        wa=source_reliability(d,['HCT116'])[0];wb=source_reliability(p,['HCT116'])[0];np.testing.assert_array_equal(wa,wb)
        self.assertEqual(wa[3],0.)
    def test_invalid_sources_abstain(self):
        rng=np.random.default_rng(52);y=rng.normal(0,.2,(5,2,160));var=np.full_like(y,.03);valid=np.ones_like(y,bool);valid[4]=False
        a=pool(y,var,valid,np.ones(5),np.array([0,1]),True);y[4]=np.nan;var[4]=np.nan
        b=pool(y,var,valid,np.ones(5),np.array([0,1]),True);np.testing.assert_array_equal(a,b)
        np.testing.assert_array_equal(a[np.arange(2),[0,1]],0.)
    def test_weights_are_bounded_and_convex(self):
        for w in GRID:self.assertGreaterEqual(w[0],.5);self.assertAlmostEqual(sum(w),1.);self.assertTrue(all(x>=0 for x in w))
    def test_quantile_transport_replays_ties_and_masks(self):
        rng=np.random.default_rng(58);x=sparse.csr_matrix(rng.poisson(.8,(30,160)));libs=np.asarray(x.sum(1)).ravel();p=moments(x)[1]
        delta=rng.normal(0,.1,(2,160,5));mask=np.arange(160)%2==0
        a=transport_initial(x,libs,p,p,delta,mask,11);b=transport_initial(x,libs,p,p,delta,mask,11)
        np.testing.assert_array_equal(a,b);self.assertTrue(np.isfinite(a).all());self.assertTrue((a>=0).all())
        np.testing.assert_allclose(a[:,~mask],x.toarray()[:,~mask],atol=1e-12)
    def test_integer_transport_preserves_supplied_libraries(self):
        rng=np.random.default_rng(59);x=sparse.csr_matrix(rng.poisson(8,(30,160)).astype(np.int32));prepared=moments(x)
        from evaluate_atlas_transfer import load_atlas_model
        atlas,_=load_atlas_model()
        count,diag=emit(x,np.arange(30),np.zeros(160),None,61,'quantile25',np.zeros((2,160,5)),np.ones(160,bool),prepared,atlas)
        np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),prepared[0]);self.assertTrue(np.issubdtype(count.dtype,np.integer));self.assertTrue((count.data>=0).all())
if __name__=='__main__':unittest.main(verbosity=2)
