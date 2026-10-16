"""Meaningful transfer invariants and masked fitting checks."""
import unittest,numpy as np
import v32_methods as m
class Transfer(unittest.TestCase):
 def test_identity_partition(self):
  labels=['GENE'+str(i) for i in range(1000)];a={t for t in labels if m.split(t)==0};b=set(labels)-a;self.assertFalse(a&b);self.assertGreater(len(a),100)
 def test_rank_scale_invariant(self):
  x=np.array([0.,1.,2.,2.,8.]);np.testing.assert_array_equal(m.control_features(x,np.arange(5)),m.control_features(7*x,np.arange(5)))
 def test_own_mask(self):
  rng=np.random.default_rng(32);x=rng.normal(size=(12,4));g=rng.normal(size=(9,3));y=rng.normal(size=(12,9));own=np.arange(12)%9;w=np.ones(12);a=m.old.masked_bilinear(x,g,y,own,w,2.);y[np.arange(12),own]+=10000;b=m.old.masked_bilinear(x,g,y,own,w,2.);np.testing.assert_allclose(a,b,atol=1e-8)
 def test_missing_target_no_tokens(self):
  # Inference loops only over observed source-target identities; absent labels
  # cannot become measured-zero response tokens.
  d=m.sources();self.assertNotIn('__unmeasured__',set(d['RPE1_full']['targets']))
if __name__=='__main__':unittest.main()