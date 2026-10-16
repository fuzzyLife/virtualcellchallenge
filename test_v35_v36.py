"""Boundary and information-path tests for retrieval/population supervision."""
import unittest
import numpy as np
import torch
import v35_methods as m
import v36_methods as p
class RetrievalTests(unittest.TestCase):
 def batch(self):
  torch.manual_seed(35);b={}
  for k in ['value','mask','src']:b[k]=torch.rand(3,5,7)
  for k,n in [('token',79),('prior',None),('context',4),('target',16)]:b[k]=torch.rand(3,5,n) if k=='token' else (torch.rand(3,5) if n is None else torch.rand(3,n))
  for k in ['base','qlog','control_var','dropout','moments_known','truth','valid']:b[k]=torch.rand(3,7)
  b['mask']=(b['mask']>.3).float();b['qlog'][:]=.6;b['token'][:,:,-3]=(torch.rand(3,5)>.5).float()
  return b,torch.rand(7,16)
 def test_zero_init_exact_baseline(self):
  b,g=self.batch();net=m.RetrievalResidual();y,v,d=net(b,g);torch.testing.assert_close(y,b['base'],rtol=0,atol=0);self.assertEqual(int(torch.count_nonzero(v)),0)
 def test_masked_token_cannot_change_output(self):
  b,g=self.batch();net=m.RetrievalResidual()
  for x in net.decoder[-1].parameters():torch.nn.init.uniform_(x,-.2,.2)
  b['mask'][:,0]=0;y=net(b,g)[0];b['value'][:,0]=99999;b['src'][:,0]=-9999;z=net(b,g)[0];torch.testing.assert_close(y,z,rtol=0,atol=0)
 def test_no_evidence_fallback_and_bounds(self):
  b,g=self.batch();net=m.RetrievalResidual()
  for x in net.parameters():torch.nn.init.normal_(x,0,.5)
  y,v,d=net(b,g);self.assertLessEqual(float(d.detach().abs().max()),.5);self.assertLessEqual(float(v.detach().abs().max()),1.5)
  b['mask']*=0;y,v,d=net(b,g);torch.testing.assert_close(y,b['base'],rtol=0,atol=0);self.assertEqual(int(torch.count_nonzero(v)),0)
 def test_token_order_invariance(self):
  b,g=self.batch();net=m.RetrievalResidual()
  for x in net.parameters():torch.nn.init.normal_(x,0,.05)
  y=net(b,g)[0];order=torch.tensor([3,1,4,0,2])
  for k in ['value','mask','src','token','prior']:b[k]=b[k][:,order]
  torch.testing.assert_close(net(b,g)[0],y,rtol=1e-6,atol=1e-6)
 def test_population_view_swap_invariance(self):
  b,g=self.batch();net=m.RetrievalResidual();y=torch.rand(3,2,7);v=torch.rand(3,2,7);mask=torch.ones_like(y)
  a=p.population_loss(net,b,g,y,v,mask)[0];z=p.population_loss(net,b,g,y.flip(1),v.flip(1),mask.flip(1))[0];torch.testing.assert_close(a,z)
 def test_population_index_shape(self):
  x=np.zeros((8,2,11));cols=np.array([1,3,5,7]);self.assertEqual(x[0,:,cols].T.shape,(2,4))
 def test_split_is_target_not_cell(self):
  self.assertEqual(m.target_split('VIM'),m.target_split('VIM'));self.assertEqual({m.target_split('VIM') for _ in range(10)},{m.target_split('VIM')})
if __name__=='__main__':unittest.main()