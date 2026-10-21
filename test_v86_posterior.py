"""Nontrivial Gaussian conditioning, missing-data and synthetic transfer checks."""
from unittest.mock import patch
import numpy as np
import v86_joint_posterior as m
from v57_pie import write,sha

def main():
    y=np.array([[.3,-.2]]);var=np.array([[.04,.09]]);ok=np.ones_like(y,bool);sim=np.array([[1.,.7],[.7,1.]])
    with patch.object(m,'COMPONENTS',[(.4,999.)]):
        ll,mean,variance=m.likelihood(y,var,ok,sim,np.array([.8,.9]));cov=.4**2*(.8*np.ones((3,3))+.2*np.eye(3));c=cov[:2,:2]+np.diag(var[0]);cross=cov[2,:2];expected=cross@np.linalg.solve(c,y[0]);expectedvar=cov[2,2]-cross@np.linalg.solve(c,cross)
        np.testing.assert_allclose(mean[0,0],expected,atol=1e-12);np.testing.assert_allclose(variance[0,0],expectedvar,atol=1e-12)
        poison=np.array([[.3,1e30]]);vpoison=np.array([[.04,np.nan]]);mask=np.array([[True,False]]);a=m.likelihood(poison,vpoison,mask,sim,np.array([.8,.9]));b=m.likelihood(np.array([[.3,-1e30]]),vpoison,mask,sim,np.array([.8,.9]))
        for x,z in zip(a,b):np.testing.assert_array_equal(x,z)
        _,zero,_=m.likelihood(np.zeros_like(y),var,ok,sim,np.array([.8,.9]));np.testing.assert_array_equal(zero,np.zeros_like(zero))
    rng=np.random.default_rng(86086);cov=.4**2*(.8*np.ones((3,3))+.2*np.eye(3));truth=rng.multivariate_normal(np.zeros(3),cov,2000);observed=truth[:,:2]+rng.normal(0,.1,(2000,2));vv=np.full((2000,2),.01);mask=np.ones_like(observed,bool)
    ll,_,_=m.likelihood(observed[:1000],vv[:1000],mask[:1000],sim);w,history=m.em(ll);ll,means,variances=m.likelihood(observed[1000:],vv[1000:],mask[1000:],sim,np.array([.7,.7]));prediction,uncertainty=m.mixture(ll,means,variances,w);ratio=float(np.mean((prediction-truth[1000:,2])**2)/np.mean((observed[1000:,0]-truth[1000:,2])**2));assert ratio<1 and (uncertainty>=0).all();np.testing.assert_allclose(w.sum(),1,atol=1e-12)
    write(m.A/'synthetic_tests.json',dict(passed=True,Gaussian_conditioning_matches_direct_solve=True,unmeasured_source_poison_has_no_effect=True,zero_observed_response_mean_exact=True,synthetic_unseen_context_MSE_ratio_to_copy=ratio,weights_sum_to_one=True,code_sha256=sha(__file__),model_code_sha256=sha(m.ROOT/'v86_joint_posterior.py')))
    print('V86 synthetic tests passed; held-context MSE ratio',ratio,flush=True)
if __name__=='__main__':main()
