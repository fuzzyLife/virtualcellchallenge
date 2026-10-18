"""Magnitude calibration preserves source identities and cannot invent signs/zeros."""
import numpy as np
from v64_response_magnitude import calibrate
rng=np.random.default_rng(64);genes=np.array([str(i) for i in range(100)]);targets=genes[:10];base=rng.normal(0,.2,(10,100)).astype(np.float32);base[:,20:24]=0;base[np.arange(10),np.arange(10)]=0
out,scale=calibrate(base,base*.7,genes,targets);np.testing.assert_allclose(scale,.7,atol=2e-7);np.testing.assert_allclose(out,base*.7,atol=2e-7);assert np.all(out[base==0]==0) and np.all(out*base>=0)
z,s=calibrate(np.zeros_like(base),np.ones_like(base),genes,targets);assert not z.any() and (s==1).all()
print('V64 scale recovery, sign/zero preservation and empty-response fallback passed.')
