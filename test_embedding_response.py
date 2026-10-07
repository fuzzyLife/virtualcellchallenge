"""Focused tests for leakage boundaries and missing-label regression behavior."""
import unittest
import numpy as np
from embedding_response_methods import kernels, predict_ridge, rank_context_features


class ResponseModelTests(unittest.TestCase):
    def test_missing_response_is_not_an_imputed_zero(self):
        k=np.eye(6);cross=np.eye(6)
        y=np.ones((6,2));v=np.ones((6,2),bool);v[0,1]=False
        a,_=predict_ridge(k,cross,y,v,np.array(['A']*3+['B']*3),ridge=1)
        y[0,1]=1e9
        b,_=predict_ridge(k,cross,y,v,np.array(['A']*3+['B']*3),ridge=1)
        np.testing.assert_array_equal(a,b)
        self.assertEqual(a[0,1],0)

    def test_pairwise_rank_retains_observed_zero(self):
        profile=np.arange(110,dtype=float)
        gap=profile.copy();gap[-5:]=np.nan
        features,n=rank_context_features(np.stack([profile,gap]),[0])
        self.assertEqual(n[1,0],105)
        self.assertAlmostEqual(features[1,0],1)
        self.assertEqual(n[0,0],110)

    def test_training_kernel_does_not_depend_on_query_control(self):
        profiles=np.array([np.arange(120),np.arange(120)[::-1],np.arange(120)],dtype=float)
        vectors=np.array([[1.,2.],[3.,4.],[5.,6.],[7.,8.],[9.,10.]])
        emb={'pca':(vectors,np.ones(5,bool)),'tahoe':(vectors,np.ones(5,bool))}
        args=(np.arange(5),np.array([0,0,0,1,1]),np.array([0]),np.array([2]))
        a,c=kernels(*args,profiles,emb,'tahoe_rank')
        profiles[2]=np.sin(np.arange(120))
        b,d=kernels(*args,profiles,emb,'tahoe_rank')
        np.testing.assert_array_equal(a,b)
        self.assertFalse(np.allclose(c,d))
        self.assertGreaterEqual(np.linalg.eigvalsh(a).min(),-1e-10)

    def test_all_unavailable_gene_is_marked_unsupported(self):
        k=np.ones((6,6));v=np.ones((6,2),bool);v[:,1]=False
        pred,mask=predict_ridge(k,np.ones((2,6)),np.ones((6,2)),v,np.array(['A']*6))
        self.assertFalse(mask[1]);np.testing.assert_array_equal(pred[:,1],0)


if __name__=='__main__':unittest.main()
