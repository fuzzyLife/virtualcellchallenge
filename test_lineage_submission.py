"""Focused tests for V4 eligible-gene correction, leakage and response policies."""
import unittest
import numpy as np
from scipy.stats import false_discovery_control
from fit_lineage_submission import bh_valid, response_prior, choose, assemble_effects

class LineageTests(unittest.TestCase):
    def test_bh_excludes_invalid_genes(self):
        p=np.array([[.001,.03,.9,.001],[.1,.001,.02,.7]])
        valid=np.array([[True,True,False,False],[True,False,True,True]])
        q=bh_valid(p,valid)
        for i in range(2):np.testing.assert_allclose(q[i,valid[i]],false_discovery_control(p[i,valid[i]]),atol=1e-14)
        np.testing.assert_array_equal(q[~valid],np.ones((~valid).sum()))

    def test_prior_never_reads_held_targets(self):
        values=np.array([[1.,2.],[3.,4.],[999.,999.]])
        errors=np.ones_like(values);valid=np.ones_like(values,dtype=bool);conditions=np.repeat('Rest',3)
        train=np.array([True,True,False])
        actual=response_prior(values,errors,valid,conditions,train)
        values[2]=[-999.,-999.]
        after=response_prior(values,errors,valid,conditions,train)
        np.testing.assert_array_equal(actual['Rest'],after['Rest'])
        np.testing.assert_array_equal(actual['Stim8hr'],[0,0])

    def test_selection_guardrail_and_ties(self):
        records=[{'setting':{'specific_gain':1,'prior_gain':0},'relative_zero_loss':1,'selection_score':.1},
          {'setting':{'specific_gain':1.5,'prior_gain':0},'relative_zero_loss':1.03,'selection_score':.9}]
        self.assertEqual(choose(records),records[0])
        transfer={'setting':{'cd4_coefficient':0,'k562_coefficient':1,'prior_gain':0},'relative_zero_loss':1,'selection_score':.3}
        self.assertEqual(choose([transfer]),transfer)

    def test_lineage_support_direct_h1_and_expression_gates(self):
        genes=np.array(['T','G','H']);targets=np.array(['T'])
        cd4=np.ones((3,1,3),dtype=np.float32);specific=cd4.copy()
        cm=np.tile([True,True,False],(3,1,1));k=np.array([[2.,2.,2.]])
        km=np.ones((1,3),bool);h=np.array([[3.,3.,3.]])
        hp=np.array([[.2,.2,.2]]);hm=np.array([[False,True,False]])
        means=np.array([[10,10,10],[10,10,10],[2.5,2.5,0]])
        effect,_=assemble_effects(cd4,specific,cm,k,km,h,hp,hm,np.array([True]),
             {'cd4_coefficient':.5,'k562_coefficient':.5},means,genes,targets)
        np.testing.assert_allclose(effect[0,0],[0,1,2])
        np.testing.assert_allclose(effect[1,0],[0,3.,1.7])
        np.testing.assert_allclose(effect[2,0],[0,.8,0])
        np.testing.assert_array_equal(cd4,np.ones((3,1,3)))

if __name__=='__main__':unittest.main()
