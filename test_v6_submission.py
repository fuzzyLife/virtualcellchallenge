"""Focused V6 safeguards: ranked DE, uncertainty correction, response and count invariants."""
import unittest
import numpy as np
from scipy import sparse,stats
from threadpoolctl import threadpool_limits
from v6_methods import moments,emit,de,direction_metrics,profile_jackknife,CANDIDATES
from fit_v6_submission import assemble,select_bounded
from evaluate_atlas_transfer import load_atlas_model
import pandas as pd


class V6Tests(unittest.TestCase):
    def test_rank_test_matches_scipy_with_ties_and_control_filter(self):
        rng=np.random.default_rng(23)
        c=rng.integers(0,30,(37,7));x=rng.integers(0,40,(31,7));c[:,0]=0;x[:,0]=100
        actual=de(sparse.csr_matrix(x),sparse.csr_matrix(c),np.array(list('ABCDEFG')))
        cn=c*1e6/c.sum(axis=1)[:,None];xn=x*1e6/x.sum(axis=1)[:,None]
        eligible=cn.mean(axis=0)>5
        expected=stats.mannwhitneyu(xn[:,eligible],cn[:,eligible],axis=0,method='asymptotic',use_continuity=True).pvalue
        np.testing.assert_allclose(actual['p'][eligible],expected,atol=1e-14)
        np.testing.assert_allclose(actual['q'][eligible],stats.false_discovery_control(expected,method='bh'),atol=1e-14)
        self.assertFalse(actual['eligible'][0]);self.assertEqual(actual['q'][0],1.)

    def test_fidelity_counts_budget_not_only_sign_precision(self):
        truth={'eligible':np.ones(3,bool),'q':np.repeat(.01,3),'lfc':np.ones(3)}
        pred={'q':np.array([.01,.1,.001]),'p':np.array([.001,.05,.0001]),'lfc':np.array([1.,-1.,1.])}
        result=direction_metrics(pred,truth,np.array(['A','B','C']),2)
        self.assertEqual(result['fidelity'],.5);self.assertEqual(result['reach'],.5);self.assertEqual(result['jaccard'],.5)
        self.assertEqual(result['predicted_calls'],1);self.assertEqual(result['reference_calls'],2)

    def test_jackknife_matches_explicit_leave_one_out(self):
        x=np.array([[3,7,1],[4,3,8],[1,1,6]],dtype=int)
        profile,jack=profile_jackknife(x)
        values=np.array([np.log1p(50000*np.delete(x,i,axis=0).sum(axis=0)/np.delete(x,i,axis=0).sum()) for i in range(3)])
        np.testing.assert_allclose(profile,np.log1p(50000*x.sum(axis=0)/x.sum()))
        np.testing.assert_allclose(jack,2/3*((values-values.mean(axis=0))**2).sum(axis=0))

    def test_response_updates_cannot_flip_sign_or_create_unsupported_effect(self):
        old=np.array([[[1.,-1.,0.,2.]],[[1.,-1.,0.,2.]],[[1.,-1.,0.,2.]]],dtype=np.float32)
        new=np.array([[[-2.,4.,1.,1.]],[[2.,2.,1.,1.]],[[2.,2.,1.,1.]]],dtype=np.float32)
        result=assemble(old,new,{'blend':.25})
        np.testing.assert_array_equal(result[0],old[0])
        np.testing.assert_array_equal(result[1,0],[1.25,-1.,0.,1.75])
        self.assertFalse((result*old<0).any());self.assertEqual(result[old==0].sum(),0.)
        with self.assertRaises(ValueError):assemble(old,new,{'blend':1.})

    def test_selection_excludes_full_replacement_and_fidelity_regression(self):
        table=pd.DataFrame([{'candidate':'v4','fidelity':.5,'jaccard':.1,'expression_ratio':1.},
            {'candidate':'v5_dual','fidelity':.9,'jaccard':.9,'expression_ratio':.1},
            {'candidate':'v4_cal10','fidelity':.4,'jaccard':.2,'expression_ratio':.5},
            {'candidate':'blend10_v4','fidelity':.51,'jaccard':.12,'expression_ratio':.9}])
        self.assertEqual(select_bounded(table),'blend10_v4')
        table.loc[3,'fidelity']=.3
        self.assertEqual(select_bounded(table),'v4')

    def test_cached_baselines_do_not_change_stochastic_emission(self):
        rng=np.random.default_rng(29);x=sparse.csr_matrix(rng.integers(10,1000,(24,12)),dtype=np.int32)
        depth,p,b,_=moments(x);prepared=(depth,p,b,np.asarray(x.mean(axis=0)).ravel())
        effect=np.linspace(-.2,.2,12);rows=np.arange(20)
        a,_=emit(x,rows,effect,1,99,CANDIDATES['v4'],None)
        b,_=emit(x,rows,effect,1,99,CANDIDATES['v4'],None,prepared=prepared)
        self.assertEqual((a-b).nnz,0)

    def test_partial_calibration_keeps_stochastic_library_sizes(self):
        rng=np.random.default_rng(31);x=sparse.csr_matrix(rng.integers(10,1000,(30,12)),dtype=np.int32);rows=np.arange(20)
        atlas,_=load_atlas_model();effect=np.linspace(-.15,.15,12)
        base,_=emit(x,rows,effect,1,17,CANDIDATES['v4'],atlas)
        corrected,diag=emit(x,rows,effect,1,17,CANDIDATES['v4_cal10'],atlas)
        replay,_=emit(x,rows,effect,1,17,CANDIDATES['v4_cal10'],atlas)
        np.testing.assert_array_equal(np.asarray(base.sum(axis=1)).ravel(),np.asarray(corrected.sum(axis=1)).ravel())
        self.assertEqual((replay-corrected).nnz,0);self.assertTrue(np.issubdtype(corrected.dtype,np.integer))
        self.assertLess(diag['requested_bulk_L1'],.002)


if __name__=='__main__':
    with threadpool_limits(limits=2):unittest.main()
