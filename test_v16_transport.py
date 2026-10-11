"""Scientific invariants for V16; no hidden competition labels."""
import unittest
import numpy as np
from scipy import sparse
from v16_methods import source_arrays,fit_calibration,predict,balanced_counts,residual_initial,emit
from v6_methods import moments
import test_v14_context as prior_tests
class TransportTests(unittest.TestCase):
    def fixture(self):
        rng=np.random.default_rng(1601);x=rng.negative_binomial(2,.2,(36,90)).astype(np.int32);x[rng.random(x.shape)<.5]=0
        return sparse.csr_matrix(x)
    def test_zero_response_is_bitwise_empirical_control(self):
        x=self.fixture();count,_=emit(x,np.arange(36),np.zeros(90),None,9,'balanced',moments(x))
        np.testing.assert_array_equal(count.toarray(),x.toarray())
    def test_transport_preserves_libraries_and_zeros(self):
        x=self.fixture();f=np.exp2(np.linspace(-2,2,90));a,e=balanced_counts(x,f,161);b,_=balanced_counts(x,f,161)
        np.testing.assert_array_equal(a,b);np.testing.assert_array_equal(a.sum(1),np.asarray(x.sum(1)).ravel());self.assertTrue((a>=0).all());self.assertTrue(np.issubdtype(a.dtype,np.integer));self.assertTrue((a[x.toarray()==0]==0).all())
    def test_expected_transport_mean_matches_redistribution(self):
        x=np.tile(np.array([40,30,20,10,0]),(6000,1));c,e=balanced_counts(sparse.csr_matrix(x),np.array([.4,.8,1.5,2.,3.]),17)
        np.testing.assert_allclose(c.mean(0),e.mean(0),atol=.25,rtol=0);np.testing.assert_array_equal(c.sum(1),100)
    def test_own_target_repression_moves_transcripts(self):
        x=sparse.csr_matrix(np.tile(np.array([100,40,60]),(3000,1)));c,_=emit(x,np.arange(3000),np.zeros(3),0,7,'balanced',moments(x))
        self.assertLess(c[:,0].mean(),25);np.testing.assert_array_equal(np.asarray(c.sum(1)).ravel(),200)
    def test_residual_preserves_poisson_coefficient(self):
        mu=np.arange(1.,31.).reshape(5,6);emp=mu[::-1].copy();noise=mu+2
        for alpha in [.1,.25]:np.testing.assert_allclose(residual_initial(mu,noise+3,emp,alpha)-residual_initial(mu,noise,emp,alpha),.65*3,atol=1e-14)
    def test_outer_held_labels_do_not_change_calibration(self):
        d=prior_tests.ContextTests().fixture();a=fit_calibration(d,['HCT116'],provider=source_arrays);p={k:v.copy() for k,v in d.items()};rows=p['row_families']=='HCT116';p['truth'][rows]=1e6;p['errors'][rows]=1e-6;p['valid'][rows]=False;p['knockdown'][rows]=.001
        b=fit_calibration(p,['HCT116'],provider=source_arrays);np.testing.assert_array_equal(a['gains'],b['gains']);np.testing.assert_array_equal(a['numerator'],b['numerator']);self.assertNotIn('HCT116',a['training_families'])
    def test_reserved_labels_do_not_fit_gain(self):
        d=prior_tests.ContextTests().fixture();a=fit_calibration(d,['HCT116'],['T1'],source_arrays);p={k:v.copy() for k,v in d.items()};r=p['row_targets']==1;p['truth'][r]=1e6;p['errors'][r]=1e-6
        b=fit_calibration(p,['HCT116'],['T1'],source_arrays);np.testing.assert_array_equal(a['gains'],b['gains']);np.testing.assert_array_equal(a['denominator'],b['denominator'])
    def test_calibration_is_bounded_and_never_reverses(self):
        d=prior_tests.ContextTests().fixture();bank=source_arrays(d,['HCT116']);m=fit_calibration(d,['HCT116'],provider=source_arrays);pred,f=predict(bank,m,d['profiles'][3])
        self.assertTrue(((f>=.5)&(f<=1)).all());self.assertTrue((abs(pred['calibrated'])<=abs(pred['v13'])+1e-7).all());self.assertTrue((pred['v13']*pred['calibrated']>=0).all());np.testing.assert_array_equal(pred['calibrated'][np.arange(2),bank['own']],0.)
    def test_depth_axis_mismatch_is_rejected(self):
        x=self.fixture()
        with self.assertRaises((IndexError,AssertionError)):emit(x,np.arange(36),np.zeros(90),None,9,'balanced',(np.ones(36)*100,*moments(x)[1:3]))
if __name__=='__main__':unittest.main(verbosity=2)
