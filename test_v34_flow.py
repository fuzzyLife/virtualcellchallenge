import unittest,numpy as np,torch
from v34_methods import Flow,integrate,coupled
from run_v34_experiment import rounded
class Transition(unittest.TestCase):
 def test_zero_velocity_identity(self):
  net=Flow(7,3)
  for v in net.parameters():torch.nn.init.zeros_(v)
  x=torch.randn(5,7);c=torch.randn(5,3);torch.testing.assert_close(integrate(net,x,c),x)
 def test_own_label_does_not_change_coupling(self):
  rng=np.random.default_rng(34);x=rng.normal(size=(8,10));y=rng.normal(size=(8,10));a=coupled(x,y,3);y[:,3]+=1000;b=coupled(x,y,3);np.testing.assert_array_equal(a,b);np.testing.assert_array_equal(a[:,3],x[:,3])
 def test_coupling_centers_follow_order(self):
  x=np.array([[0.,0.],[10.,0.]]);y=x[::-1].copy();a,order=coupled(x,y,-1,True);np.testing.assert_array_equal(a,y[order])
 def test_rounding_preserves_cell_libraries(self):
  q=np.random.default_rng(34).uniform(size=(4,30));libs=np.array([53,60,99,74]);x=rounded(q,libs);np.testing.assert_array_equal(x.sum(1).A1,libs);self.assertTrue((x.data>=0).all())
if __name__=='__main__':unittest.main()