"""Check compensation, bounds and zero-response identity for V45."""
import numpy as np
import v45_bounded_dispersion as m
def main():
 e=np.array([-3.,-.2,-.1,0.,.1,.2,3.])
 for arm in m.ARMS:
  w=m.residual_weight(e,arm);assert (w>=.1).all() and (w<=.2).all();np.testing.assert_allclose(w,w[::-1])
 assert m.residual_weight(e,'response_weighted20')[3]==.1
 rng=np.random.default_rng(45);x=rng.random((4,7))
 for arm in m.ARMS:np.testing.assert_allclose(m.old.initial(x,x,x,m.residual_weight(e,arm)),x,atol=1e-15)
 m.A.mkdir(parents=True,exist_ok=True);m.write(m.A/'contract_tests.json',{'passed':True,'checks':['coefficients10–20%','sign symmetry','zero-response weighted arm identical10%','compensated equal-component identity'],'script_sha256':m.b.sha256(m.ROOT/'test_v45_bounded_dispersion.py')})
 print('V45 contract tests PASS',flush=True)
if __name__=='__main__':main()

