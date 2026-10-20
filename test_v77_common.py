"""V77 common-component identities, norm-matched controls and bounded changes."""
import numpy as np
from v77_common_response import calibrate,permute_common
rng=np.random.default_rng(7777);q=rng.dirichlet(np.ones(40));q[-1]=0;q/=q.sum();own=np.arange(8);e=rng.uniform(-.4,.4,(8,40)).astype(np.float32);e[np.arange(8),own]=0;e[0]=0
for arm in ('common25','common50','common100','projection50'):
 for shuffled in (False,True):
  out=calibrate(e,q,own,arm,shuffled);assert np.isfinite(out).all() and np.max(abs(out-e))<=.375+1e-6
  np.testing.assert_array_equal(out[0],e[0]);np.testing.assert_array_equal(out[:,-1],e[:,-1]);assert (out[np.arange(8),own]==0).all()
panel=np.zeros(40,bool);panel[own]=True;c=rng.normal(size=(8,40));permuted=permute_common(c,q,panel)
np.testing.assert_allclose(np.sum(c[:,~panel]**2,axis=1),np.sum(permuted[:,~panel]**2,axis=1),rtol=1e-14)
print('V77 no-new-zero-response, absent-gene, own handling, bounds and control-norm checks passed.')
