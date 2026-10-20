"""Moment calibration with negative parent dispersion and exact output margins."""
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
from v79_count_innovations import innovation_parameters,draw_innovations,transform
from v69_count_bridge import parameters
rng=np.random.default_rng(7979);mu=np.ones((100000,1));obs=rng.binomial(2,.5,size=(len(mu),1));a0=np.array([-.5]);target=np.array([.5])
for keep in (0.,.5,.8):
 rho,add=innovation_parameters(a0,target,keep)
 variance=rho*rho*.5+rho*(1-rho)+(1-rho)+add*(1-rho)**2
 np.testing.assert_allclose(variance,1.5,rtol=0,atol=1e-12)
 y=draw_innovations(obs,mu,rho,add,rng);assert abs(y.mean()-1)<.03 and abs(y.var()-1.5)<.15
rho,add=innovation_parameters(np.array([9.]),np.array([.01]),.8);assert add[0]>=0 and rho[0]<.8
p=rng.dirichlet(np.full(96,.8));libs=rng.integers(2000,10000,48);x=sparse.csr_matrix(np.array([rng.multinomial(int(n),p) for n in libs],dtype=np.int32));q=parameters(x);q=dict(q,alpha=np.clip(q['alpha']*2,.005,20))
with threadpool_limits(limits=2):
 out,diag=transform(x,q,np.zeros(96),'turnover50',7979)
 np.testing.assert_array_equal(out.sum(0),x.sum(0));np.testing.assert_array_equal(out.sum(1),x.sum(1));assert diag['CPM_L1']<=.0005
print('V79 underdispersion/innovation moment, stochastic, adapted-retention and exact-margin checks passed.',flush=True)
