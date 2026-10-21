"""Operator arithmetic checks using independent low-rank synthetic interventions."""
import numpy as np
from threadpoolctl import threadpool_limits
from v57_pie import write,sha
import v84_intervention_operator as m
rng=np.random.default_rng(84084);latent=rng.normal(size=(100,8));loading=rng.normal(size=(8,128));x=latent@loading*.2;effects=np.array([x,1.2*x,.7*x,1.5*x]);payload=dict(effects=effects,targets=np.array([f'T{i}' for i in range(100)]));original=m.pairs;m.pairs=lambda:payload
with threadpool_limits(limits=2):
 fit=m.fit(np.arange(80));fake=m.fit(np.arange(80),True)
 for kind in m.ARMS:
  real=m.transform(x[80:],fit,kind);bad=m.transform(x[80:],fake,kind);truth=effects[1:,80:]-x[None,80:];identity=float(np.mean(truth**2));error=float(np.mean((real-truth)**2));fake_error=float(np.mean((bad-truth)**2));assert error<identity*.25 and error<fake_error,(kind,error,identity,fake_error);np.testing.assert_array_equal(m.transform(np.zeros((2,128)),fit,kind),np.zeros((3,2,128)));print('Synthetic',kind,'ratio',error/identity,flush=True)
m.pairs=original;write(m.A/'numerical_tests.json',dict(passed=True,low_rank_known_operator_recovery=True,zero_response_maps_exactly_zero=True,fitted_pair_derangement_control=True,script_sha256=sha(__file__)))
