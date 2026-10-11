"""Scientific invariants for V14 interactions; no hidden competition labels."""
import unittest
import numpy as np
from v14_context_methods import ridge_interactions,feature_fit,feature_query,fit,predict


class ContextTests(unittest.TestCase):
    def test_known_context_reversal(self):
        z=np.array([[-1.],[0.],[1.]])
        y=z[:,None,:];v=np.ones_like(y)*.01;m=np.ones_like(y,bool)
        beta,center,support=ridge_interactions(y,v,m,z,.1)
        self.assertGreater(float(beta[0,0,0]),.9)
        self.assertLess(float(np.sum(beta[0,0]*(-1-center[0,0]))),0)
        self.assertGreater(float(np.sum(beta[0,0]*(1-center[0,0]))),0)

    def test_missing_values_do_not_vote_zero(self):
        z=np.array([[-1.],[0.],[1.],[2.]])
        y=np.array([-1.,0.,1.,1e30])[:,None,None];v=np.ones_like(y)*.01;m=np.ones_like(y,bool);m[-1]=False
        a=ridge_interactions(y,v,m,z,2.)
        y[-1]=np.nan;v[-1]=np.nan
        b=ridge_interactions(y,v,m,z,2.)
        np.testing.assert_array_equal(a[0],b[0]);np.testing.assert_array_equal(a[1],b[1])
        m[0]=False;self.assertEqual(float(abs(ridge_interactions(y,v,m,z,2.)[0]).sum()),0.)

    def fixture(self):
        rng=np.random.default_rng(4401);genes=np.array(['T0','T1']+[f'G{i}' for i in range(158)]);targets=genes[:2]
        families=np.array(['K562','CD4','H1','HCT116','HEK293T'])
        rf=np.repeat(families,2);rt=np.tile([0,1],5);rc=np.repeat(np.arange(5),2)
        y=rng.normal(0,.2,(10,160));valid=np.ones_like(y,bool);valid[np.arange(10),rt]=False
        return {'genes':genes,'targets':targets,'truth':y.astype(np.float32),'errors':np.full_like(y,.15,dtype=np.float32),
                'valid':valid,'knockdown':np.full(10,.9),'profiles':rng.uniform(10,100,(5,160)),
                'row_families':rf,'row_targets':rt,'row_contexts':rc}

    def test_whole_family_labels_cannot_change_predictions(self):
        d=self.fixture();a=fit(d,'HCT116');query=d['profiles'][3].copy();pa,_=predict(a,query)
        poisoned={k:v.copy() for k,v in d.items()};rows=poisoned['row_families']=='HCT116'
        poisoned['truth'][rows]=1e5;poisoned['errors'][rows]=1e-6;poisoned['knockdown'][rows]=.01
        poisoned['valid'][rows]=False
        b=fit(poisoned,'HCT116');pb,_=predict(b,query)
        for name in pa:np.testing.assert_array_equal(pa[name],pb[name])
        self.assertNotIn('HCT116',a['training_families'])
        for name in pa:np.testing.assert_array_equal(pa[name][np.arange(2),[0,1]],0.)

    def test_poor_gene_coverage_abstains(self):
        d=self.fixture();m,z=feature_fit(d['profiles'],[0,1,2,4],d['genes'],d['targets'])
        q=d['profiles'][3].copy();q[:100]=np.nan
        z,diagnostic=feature_query(m,q);self.assertTrue(diagnostic['abstain']);np.testing.assert_array_equal(z,0.)


if __name__=='__main__':unittest.main(verbosity=2)
