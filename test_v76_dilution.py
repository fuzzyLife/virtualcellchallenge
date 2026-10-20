"""Complete-profile dilution conservation and identity checks."""
import numpy as np
from v76_response_dilution import probability,parent_prob
rng=np.random.default_rng(7676);q=rng.dirichlet(np.ones(24));q[-1]=0;q/=q.sum();e=rng.uniform(-.5,.5,(8,24));own=np.arange(8);e[np.arange(8),own]=0
p0=parent_prob(e,q,own)
for arm,gain in [('mixture075',.75),('mixture050',.5)]:
 p=probability(e,q,own,arm);np.testing.assert_allclose(p,(1-gain)*q+gain*p0,rtol=0,atol=1e-15)
 np.testing.assert_allclose(p.sum(1),1,atol=1e-15);assert (p>=0).all() and (p[:,-1]==0).all()
 sh=probability(e,q,own,arm,True);np.testing.assert_allclose(sh.sum(1),1,atol=1e-15);np.testing.assert_array_equal(sh[np.arange(8),own],p[np.arange(8),own])
for arm in ('dose030','dose055'):
 p=probability(e,q,own,arm);np.testing.assert_allclose(p.sum(1),1,atol=1e-15);assert (p>=0).all()
print('V76 exact convex mixture, closed mass, zero-gene and matched-own negative-control checks passed.')
