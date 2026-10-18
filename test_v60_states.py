"""V60 nuisance residualization: no treatment-fitting and useful-signal contracts."""
import numpy as np
from v60_depth_adjusted_states import residualize

def main():
 rng=np.random.default_rng(60);n=1000;cov=rng.normal(size=(n,2));signal=rng.normal(size=n);values=3*cov[:,0,None]+cov[:,1,None]+signal[:,None]+rng.normal(0,.1,size=(n,8));ntc=np.arange(800)
 v,r=residualize(values,ntc,cov);changed=values.copy();changed[800:]+=100;vv,rr=residualize(changed,ntc,cov)
 assert r==rr;np.testing.assert_array_equal(v[:800],vv[:800])
 assert abs(np.corrcoef(v[:800,0],cov[:800,0])[0,1])<.01
 assert np.corrcoef(v[:800,0],signal[:800])[0,1]>.95
 assert np.isfinite(v).all()
 print('V60 NTC-only fitting and nuisance/signal contracts passed',flush=True)
if __name__=='__main__':main()
