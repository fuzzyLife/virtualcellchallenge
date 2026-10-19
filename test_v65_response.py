"""Numerical and exclusion tests for available-axis response conversion."""
import numpy as np
import v65_closed_response as m

def main():
 rng=np.random.default_rng(65);f,t,g=3,30,40
 q=rng.uniform(1,100,(f,g));q[0,0]=np.nan;own=np.arange(t)
 bank=dict(profiles=q,y=rng.normal(0,.15,(f,t,g)),var=np.full((f,t,g),.05),ok=np.ones((f,t,g),bool),own=own)
 y,v,closed=m.source_arrays(bank,True)
 assert np.isfinite(y).all() and np.isfinite(v).all() and (v>=.0025).all()
 np.testing.assert_allclose(closed.sum(1),1e6,rtol=1e-14);assert closed[0,0]==0
 # Source closure removes arbitrary original library coverage scale.
 bank2=dict(bank,profiles=q*np.array([.3,2,10])[:,None])
 y2,v2,_=m.source_arrays(bank2,True)
 np.testing.assert_allclose(y,y2,atol=2e-14);np.testing.assert_allclose(v,v2,atol=2e-14)
 # No invalid response becomes evidence or contributes a changed-count normalization.
 bank3=dict(bank,y=bank['y'].copy(),ok=bank['ok'].copy());bank3['ok'][0,:,2]=False
 first=m.source_arrays(bank3,True);bank3['y'][0,:,2]=99;second=m.source_arrays(bank3,True)
 np.testing.assert_array_equal(first[0],second[0]);np.testing.assert_array_equal(first[1],second[1])
 assert not np.array_equal(m.source_arrays(bank,False)[0],y)
 print('V65 passed finite values, available-axis closure, scale invariance, invalid-response poisoning and nontrivial-composition tests.')
if __name__=='__main__':main()
