import unittest,numpy as np,torch
from v34_corrected_methods import integrate
class Constant(torch.nn.Module):
 def forward(self,x,x0,t,c):return torch.ones_like(x)
class OwnLock(unittest.TestCase):
 def test_missing_own_output_cannot_drive_other_steps(self):
  x=torch.randn(5,7);c=torch.zeros(5,2);y=integrate(Constant(),x,c,own=3);torch.testing.assert_close(y[:,3],x[:,3]);keep=np.arange(7)!=3;torch.testing.assert_close(y[:,keep],x[:,keep]+1)
if __name__=='__main__':unittest.main()