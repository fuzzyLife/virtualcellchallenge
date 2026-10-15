"""V26 focused likelihood, inference and integer-count contract tests."""
import unittest
import numpy as np,torch
from v26_network import Decoder
from v26_methods import decode,fit_probabilities,integerize,repair_cpm
class Counts(unittest.TestCase):
 def test_network_gradients_and_numpy_replay(self):
  torch.manual_seed(26);n=Decoder(12,31);x=torch.randn(7,12);s=torch.randn(7,8);c=torch.randn(7,16);base=torch.softmax(torch.randn(7,31),1);lp,kl=n(x,s,c,base,torch.ones(31));loss=-(base*lp).sum(1).mean()+.001*kl.mean();loss.backward();self.assertTrue(all(torch.isfinite(p.grad).all() for p in n.parameters()));z=torch.randn(7,16);w={k:v.detach().numpy() for k,v in n.state_dict().items()};np.testing.assert_allclose(decode({'weights':w},z.numpy(),s.numpy(),c.numpy()),n.residual(z,s,c).detach().numpy(),rtol=1e-5,atol=1e-6)
 def test_feasible_moments_and_integer_margins(self):
  rng=np.random.default_rng(26);ref=rng.poisson(20,(80,35)).astype(np.int32);libs=ref.sum(1);p=(ref/libs[:,None]).mean(0);bulk=ref.sum(0)/libs.sum();initial=rng.lognormal(0,1,(80,35));prob=fit_probabilities(initial,libs,p,bulk);x,moved=integerize(prob,libs,ref.sum(0),26);self.assertTrue(np.issubdtype(x.dtype,np.integer));np.testing.assert_array_equal(x.sum(1),libs);np.testing.assert_array_equal(x.sum(0),ref.sum(0));xx,_=integerize(prob,libs,ref.sum(0),26);np.testing.assert_array_equal(x,xx);self.assertLess(abs(prob.mean(0)-p).sum(),.01)
 def test_cpm_swaps_protect_both_margins(self):
  rng=np.random.default_rng(27);ref=rng.poisson(15,(50,25)).astype(np.int32);libs=ref.sum(1);q=rng.uniform(size=ref.shape);q/=q.sum(1,keepdims=True);x,_=integerize(q,libs,ref.sum(0),27);before=abs((x/libs[:,None]).mean(0)-(ref/libs[:,None]).mean(0)).sum();y,moved=repair_cpm(x.copy(),ref,max_fraction=.5);after=abs((y/libs[:,None]).mean(0)-(ref/libs[:,None]).mean(0)).sum();self.assertLessEqual(after,before+1e-12);np.testing.assert_array_equal(y.sum(0),ref.sum(0));np.testing.assert_array_equal(y.sum(1),libs);self.assertLessEqual(moved,int(libs.sum()*.5))
 def test_zero_gene_support_preserved(self):
  initial=np.ones((20,5));libs=np.arange(100,120);p=np.array([.1,.2,.3,.4,0.]);q=fit_probabilities(initial,libs,p,p);np.testing.assert_array_equal(q[:,-1],np.zeros(20));np.testing.assert_allclose(q.sum(1),1)
if __name__=='__main__':unittest.main()
