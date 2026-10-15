"""Independent mathematical checks for masked bilinear response fitting."""
import unittest
import numpy as np
from v30_methods import masked_bilinear
class BilinearTests(unittest.TestCase):
 def test_masked_normal_equations(self):
  r=np.random.default_rng(30);x=r.normal(size=(9,4));g=r.normal(size=(7,3));y=r.normal(size=(9,7));own=np.array([0,1,2,3,4,-1,5,6,-1]);w=r.uniform(.2,2,9);penalty=.7;rows=[];labels=[]
  for i in range(9):
   for j in range(7):
    if j==own[i]:continue
    rows.append(np.outer(x[i],g[j]).ravel()*np.sqrt(w[i]));labels.append(y[i,j]*np.sqrt(w[i]))
  xx=np.array(rows);yy=np.array(labels);expected=np.linalg.solve(xx.T@xx+penalty*np.eye(12),xx.T@yy).reshape(4,3)
  z=y.copy();ii=np.flatnonzero(own>=0);z[ii,own[ii]]=0;actual=masked_bilinear(x,g,z,own,w,penalty);np.testing.assert_allclose(actual,expected,atol=1e-7,rtol=1e-7)
 def test_unobserved_own_values_do_not_enter(self):
  r=np.random.default_rng(31);x=r.normal(size=(5,3));g=r.normal(size=(6,2));y=r.normal(size=(5,6));own=np.arange(5);y[np.arange(5),own]=0;actual=masked_bilinear(x,g,y,own,np.ones(5),1);other=y.copy();other[np.arange(5),own]=1e5;np.testing.assert_allclose(masked_bilinear(x,g,other,own,np.ones(5),1),actual,rtol=1e-6,atol=1e-6)
if __name__=='__main__': unittest.main()
