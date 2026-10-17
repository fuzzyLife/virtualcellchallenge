"""Averaged propensity equals the un-clipped mean of source confidence boosts."""
import numpy as np
import v50_response_prior as m
def main():
 base=np.array([[.2,-.2,.01,0.],[.3,-.3,.01,0.],[.4,-.4,.01,0.]],np.float32);z=np.array([[.5,-.5,np.nan,1.],[1.,-.3,2.,-1.],[-.2,.4,1.,0.]],float);sat=np.zeros(z.shape,bool);out=m.averaged(base,z,sat)
 for i in range(3):
  for j in range(2):
   other=np.arange(3)!=i;conf=np.tanh(abs(z[other,j]))*(z[other,j]*base[i,j]>0);expected=base[i,j]*(1+.25*conf.mean());np.testing.assert_allclose(out[i,j],expected,rtol=1e-6)
 np.testing.assert_array_equal(out[:,2:],base[:,2:]);assert (abs(out)<=1.25*abs(base)+1e-6).all()
 m.A.mkdir(parents=True,exist_ok=True);m.write(m.A/'contract_tests.json',{'passed':True,'checks':['leave-own-intervention-out signed propensity','mean confidence boost identity away from clipping','weak/zero unchanged','positive gains bounded1..1.25'],'script_sha256':m.b.sha256(m.ROOT/'test_v50_response_prior.py')});print('V50 response-prior contracts PASS',flush=True)
if __name__=='__main__':main()

