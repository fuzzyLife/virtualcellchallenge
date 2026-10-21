"""Independent numerical fixtures for V88 solvers and atom-aware mixture quantiles."""
import numpy as np
from scipy.optimize import minimize,brentq
from scipy.special import ndtr
from threadpoolctl import threadpool_limits
import v88_mixture_quantiles as m
from v57_pie import write

def main():
 rng=np.random.default_rng(88)
 for k in range(1,6):
  n=83;a=rng.normal(size=(k,k));cov=a@a.T/k;var=np.exp(rng.uniform(-6,3,(n,k)));y=rng.normal(size=(n,k));ll,inv,alpha=m.gaussian(y,var,cov);c=np.broadcast_to(cov,(n,k,k)).copy();ix=np.arange(k);c[:,ix,ix]+=var;expected=np.linalg.inv(c);sg,ld=np.linalg.slogdet(c);ea=np.linalg.solve(c,y[...,None])[...,0];el=-.5*(ld+np.sum(y*ea,axis=1)+k*np.log(2*np.pi));np.testing.assert_allclose(inv,expected,rtol=3e-11,atol=1e-10);np.testing.assert_allclose(alpha,ea,rtol=3e-11,atol=1e-10);np.testing.assert_allclose(ll,el,rtol=3e-11,atol=1e-10)
 ll=rng.normal(size=(500,6))*2;w,h=m.solve_weights(ll);a=np.exp(ll-ll.max(1,keepdims=True))
 fun=lambda x:-m.objective(a,x)/len(a)
 jac=lambda x:-(a.T@(1/(a@x))+.5/x)/len(a)
 fit=minimize(fun,np.full(6,1/6),jac=jac,bounds=[(1e-10,1)]*6,constraints={'type':'eq','fun':lambda x:x.sum()-1,'jac':lambda x:np.ones(6)},method='SLSQP',options={'ftol':1e-13,'maxiter':1000});assert fit.success,fit.message;np.testing.assert_allclose(w,fit.x,rtol=0,atol=2e-6)
 mu=rng.normal(size=(40,4));mu[:,0]=0;vv=rng.uniform(.001,2,(40,4));vv[:,0]=0;prob=rng.dirichlet([2,1,1,1],40)
 for q in [.1,.9]:
  result=m.quantile(mu,vv,prob,q)
  for i in range(40):
   def continuous(x):return sum(prob[i,j]*ndtr((x-mu[i,j])/np.sqrt(vv[i,j])) for j in range(1,4))
   left=continuous(0);right=left+prob[i,0]
   if left<q<=right:expected=0
   elif q<=left:expected=brentq(lambda x:continuous(x)-q,-30,0)
   else:expected=brentq(lambda x:continuous(x)+prob[i,0]-q,0,30)
   np.testing.assert_allclose(result[i],expected,rtol=0,atol=3e-8)
 # Pure-null posterior must project any nonzero contrast exactly to its atom.
 mu=np.zeros((4,1));v=np.zeros_like(mu);ll=np.zeros_like(mu);p=m.project(np.array([-2.,0.,.3,3.]),mu,v,ll,np.ones(1));np.testing.assert_array_equal(p,np.zeros(4))
 write(m.A/'synthetic_tests.json',dict(passed=True,gaussian_dimensions_checked=[1,2,3,4,5],independent_SLSQP_weight_check=True,weight_KKT_residual=h[-1]['kkt_relative_residual'],independent_Brent_quantile_checks=80,null_atom_projection_exact=True,scope='Numerical implementation tests, not biological or leaderboard validation.'));print('V88 solver and atom-aware interval fixtures passed',flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
