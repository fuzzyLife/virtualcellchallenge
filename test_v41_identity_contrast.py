"""Matched control, source-validation abstention and exact baseline tests."""
from unittest.mock import patch
import numpy as np
import v41_identity_contrast as m
def records(ce,other,before):
 return ({'sha256':'identity','before':{'cross_entropy':before},'history':[{'objective':1.,'cross_entropy':ce}]},
         {'sha256':'control','before':{'cross_entropy':before},'history':[{'objective':1.,'cross_entropy':other}]})
def main():
 for ce,other,before,expected in [(2.,2.,3.,False),(4.,5.,3.,False),(2.,3.,4.,True)]:
  a,b=records(ce,other,before)
  with patch.object(m.old,'load',side_effect=[np.ones((2,100),np.float32),np.zeros((2,100),np.float32)]),patch.object(m.old,'read',side_effect=[a,b]):
   delta,r=m.calibration('H1');assert r['use_calibrator']==expected;np.testing.assert_array_equal(delta,np.ones((2,100)))
 rng=np.random.default_rng(41);base=rng.normal(0,.3,(10,100)).astype(np.float32);base[:,0]=0;delta=rng.normal(0,2,(2,100)).astype(np.float32)
 np.testing.assert_array_equal(m.predict(delta,base,False),base);np.testing.assert_array_equal(m.predict(np.zeros_like(delta),base,True),base)
 out=m.predict(delta,base,True);assert (out*base>=0).all();assert (out[:,0]==0).all();ratio=out[:,1:]/base[:,1:];assert ratio.min()>=.75-1e-6 and ratio.max()<=1.25+1e-6
 m.write(m.A/'contract_tests.json',{'passed':True,'matched_or_worse_source_learner_abstains':True,'source_learner_must_beat_baseline_and_control':True,'abstention_and_identical_models_give_exact_baseline':True,'positive_gain_and_own_zero_contract':True,'script_sha256':m.b.sha256(m.ROOT/'test_v41_identity_contrast.py')});print('V41 matched-control and abstention contracts PASS',flush=True)
if __name__=='__main__':main()

