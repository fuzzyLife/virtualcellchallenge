"""Test exact closure, integer-emitter own convention and margin feasibility."""
import numpy as np
from v63_identity_guard import guard,probabilities
rng=np.random.default_rng(63);base=rng.normal(0,.2,(12,100));new=rng.normal(0,1.,(12,100));own=np.arange(12);base[np.arange(12),own]=0;new[np.arange(12),own]=0;bulk=rng.uniform(.1,1,100);bulk/=bulk.sum();panel=np.arange(100)<12
out,rec=guard(base,new,bulk,own,panel);assert rec['all_target_margin_constraints_passed'] and np.isfinite(out).all();assert (out[np.arange(12),own]==0).all()
p0=probabilities(base,bulk,own);pn=probabilities(new,bulk,own);actual=probabilities(out,bulk,own);alpha=np.array(rec['fractions']);np.testing.assert_allclose(actual,p0+alpha[:,None]*(pn-p0),atol=2e-9,rtol=2e-7)
same,r=guard(base,base,bulk,own,panel);np.testing.assert_allclose(same,base,atol=2e-7);assert r['full_correction_targets']==12
print('V63 exact probability closure, own knockdown convention and all target-margin constraints passed.')
