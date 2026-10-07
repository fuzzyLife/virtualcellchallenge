"""Scientific invariance checks for V7 source preparation and response pooling."""
import unittest
import numpy as np,pandas as pd
from scipy import sparse
from prepare_v7_sources import balanced_sample
from v7_response_methods import combine,adaptive_shrink
from fit_v7_response import summarize_raw
class ResponseTests(unittest.TestCase):
    def test_unmeasured_source_abstains_but_measured_zero_votes(self):
        y=np.array([[[1.,1.]],[[0.,0.]]]);v=np.ones_like(y)*.01
        mask=np.array([[[True,True]],[[False,True]]]);kd=np.ones((2,1))*.9
        out,noise,ok=combine((y,v,mask,kd),np.ones(2))
        self.assertAlmostEqual(out[0,0],1.);self.assertAlmostEqual(out[0,1],.5)
        self.assertTrue(ok.all());self.assertGreater(noise[0,1],noise[0,0])
    def test_all_missing_is_unsupported(self):
        y=np.ones((2,1,3));v=np.ones_like(y);mask=np.zeros_like(y,dtype=bool)
        out,noise,ok=combine((y,v,mask,np.ones((2,1))),np.ones(2))
        self.assertFalse(ok.any());self.assertTrue((out==0).all())
    def test_adaptive_preserves_sign_and_reduces_noise(self):
        values=np.tile(np.r_[np.zeros(100),np.ones(30),-np.ones(30)],(1,1))
        var=np.ones_like(values)*.01;valid=np.ones_like(values,bool)
        result,mix=adaptive_shrink(values,var,valid)
        self.assertTrue((np.sign(result)==np.sign(values)).all())
        self.assertTrue((abs(result)<=abs(values)+1e-7).all())
        self.assertGreater(result[0,110],.7)
        np.testing.assert_allclose(mix.sum(axis=1),1)
    def test_batch_balancing_and_no_duplicate_cells(self):
        frame=pd.DataFrame({'sample':np.repeat(['a','b'],[50,5]),'cell':np.arange(55)})
        sample=balanced_sample(frame,8,7)
        self.assertEqual(sample.groupby('sample').size().to_dict(),{'a':4,'b':4})
        self.assertTrue(sample.cell.is_unique)
        pd.testing.assert_frame_equal(sample,balanced_sample(frame,8,7))
    def test_matched_batches_do_not_create_neutral_response(self):
        genes=np.array(['target','g1','g2']);targets=np.array(['target'])
        values=[];meta=[]
        # Treated cells occur only in a; controls also occur in b, with different composition.
        for batch,label,n,count in [('a','target',40,[10,80,10]),('a','Non-Targeting',40,[10,80,10]),('b','Non-Targeting',40,[10,10,80])]:
            for i in range(n):values.append(count);meta.append({'sample':batch,'gene_target':label,'partition':0,'total_counts':100})
        result=summarize_raw(sparse.csr_matrix(values),pd.DataFrame(meta),genes,targets,np.ones(3,bool),2)
        np.testing.assert_allclose(result['truth'],0,atol=1e-7)
        self.assertFalse(result['valid'][0,0]);self.assertTrue(result['valid'][0,1:].all())
    def test_reduced_noise_preserves_count_expectation(self):
        from v7_count_methods import emit
        x=sparse.csr_matrix(np.tile([100,20,5],(400,1)),dtype=np.int32)
        effect=np.log2([.6,1.5,1.]);draws={key:[] for key in ['v4','low_noise']}
        for seed in range(32):
            for key in draws:draws[key].append(emit(x,np.arange(400),effect,None,seed,key).toarray())
        original=np.concatenate(draws['v4']);reduced=np.concatenate(draws['low_noise'])
        np.testing.assert_allclose(reduced.mean(axis=0),[60,30,5],atol=.15)
        np.testing.assert_allclose(original.mean(axis=0),[60,30,5],atol=.15)
        ratio=reduced[:,:2].var(axis=0)/original[:,:2].var(axis=0)
        self.assertTrue(((ratio>.45)&(ratio<.65)).all())
        self.assertTrue((reduced[:,2]==5).all())
    def test_neutral_emitters_preserve_template_identity(self):
        from v7_count_methods import emit
        x=sparse.csr_matrix([[0,3,9],[2,1,4]],dtype=np.int32)
        for generator in ['v4','low_noise']:
            result=emit(x,np.arange(2),np.zeros(3),None,17,generator)
            self.assertEqual((result!=x).nnz,0)
            self.assertTrue(np.issubdtype(result.dtype,np.integer))
    def test_ridge_mask_does_not_impute_missing_labels(self):
        from train_v7_ridge import rectangular_mask
        data={'genes':np.array(['target','g1','g2']),'targets':np.array(['target']),
            'row_contexts':np.array([0,0,1]),'valid':np.array([[False,True,True],[True,True,False],[True,False,True]])}
        mask=rectangular_mask(data)
        self.assertFalse((mask&~data['valid']).any())
        self.assertFalse(mask[:,0].any());self.assertTrue(mask[:2,1].all());self.assertFalse(mask[:2,2].any())
        self.assertTrue(mask[2,2])
    def test_population_moments_match_treated_batches(self):
        from evaluate_v7_response import population_statistics,matched_population
        x=sparse.csr_matrix([[10,90],[20,80],[90,10],[80,20],[1,99],[2,98]],dtype=np.int32)
        meta=pd.DataFrame({'sample':['a','a','b','b','a','a']})
        population=population_statistics(x,meta,np.arange(4));depth,p,b,raw=matched_population(population,meta,np.array([4,5]))
        np.testing.assert_allclose(raw,[15,85]);np.testing.assert_allclose(p,[.15,.85]);np.testing.assert_allclose(b,p)
        self.assertEqual(depth[0],100)
if __name__=='__main__':unittest.main(verbosity=2)

