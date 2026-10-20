"""Assay-gain calibration identities and correction bounds."""
import numpy as np
from v74_assay_gain import calibrate
rng=np.random.default_rng(7474);q=rng.dirichlet(np.ones(24));q[-1]=0;q/=q.sum();e=rng.uniform(-.3,.3,(8,24)).astype(np.float32);own=np.arange(8);e[np.arange(8),own]=0;e[0]=0
np.testing.assert_array_equal(calibrate(e,q,own,1),e)
for gain in (1.15,1.3,1.5):
 result=calibrate(e,q,own,gain)
 np.testing.assert_array_equal(result[0],e[0]);np.testing.assert_array_equal(result[:,-1],e[:,-1]);assert (result[np.arange(8),own]==0).all()
 assert np.max(abs(result-e))<=.375+1e-6 and np.isfinite(result).all()
q[own[-1]]=0;q/=q.sum();assert np.isfinite(calibrate(e,q,own,1.3)).all()
print('V74 gain-one/zero-response identities, absent-gene handling and bounds passed.')
