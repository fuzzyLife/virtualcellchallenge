import unittest,numpy as np
from run_v33_experiment import latent_labels,cap_rows
class Labels(unittest.TestCase):
 def test_missing_own_least_squares(self):
  rng=np.random.default_rng(33);u,_=np.linalg.qr(rng.normal(size=(30,8)));y=rng.normal(size=(4,30));own=np.array([0,7,12,-1]);y[np.arange(3),own[:3]]=0;g={'basis':u,'scale':np.ones(30)};a=latent_labels(g,y,own)
  for i,o in enumerate(own):
   keep=np.arange(30)!=o;b=np.linalg.lstsq(u[keep],y[i,keep],rcond=None)[0];np.testing.assert_allclose(a[i],b,atol=1e-10)
 def test_cap_excludes_validation(self):
  d={'validation':np.array([False]*600+[True]*100),'routes':np.array(['a']*300+['b']*300+['c']*100)};ids=cap_rows(d);self.assertEqual(len(ids),512);self.assertFalse(d['validation'][ids].any());self.assertEqual(len(set(ids)),512)
if __name__=='__main__':unittest.main()