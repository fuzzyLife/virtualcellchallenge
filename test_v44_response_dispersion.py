"""Conservation and bound checks for residual weighting."""
import numpy as np
import v44_response_dispersion as m
def main():
 e=np.array([-3.,-.1,0.,.1,3.]);rng=np.random.default_rng(44);expected=rng.uniform(0,20,(8,5));poisson=rng.poisson(expected);emp=rng.uniform(0,20,(8,5))
 for a in m.ARMS:
  w=m.residual_weight(e,a);assert (w>=.1).all() and (w<=.35).all();q=m.initial(expected,poisson,emp,w);assert (q>=0).all()
  # With equal components, weighting cannot alter the profile.
  np.testing.assert_allclose(m.initial(expected,expected,expected,w),expected,atol=1e-12)
 np.testing.assert_allclose(m.initial(expected,poisson,emp,np.full(5,.1)),.25*expected+.65*poisson+.1*emp,atol=1e-12)
 w=m.residual_weight(e,'adaptive');assert w[2]>w[1]>w[0] and w[0]==w[-1]
 m.A.mkdir(parents=True,exist_ok=True);m.write(m.A/'contract_tests.json',{'passed':True,'tests':['nonnegative residual coefficients','convex equal-component identity','alpha.1 reproduces V24 initial algebra','adaptive weighting symmetric and stronger for weak effects'],'script_sha256':m.b.sha256(m.ROOT/'test_v44_response_dispersion.py')});print('V44 contract tests PASS',flush=True)
if __name__=='__main__':main()

