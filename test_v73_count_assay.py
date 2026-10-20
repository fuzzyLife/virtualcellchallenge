"""Null-centered likelihood, known gain, sparse-evidence and assay-inversion tests."""
import numpy as np
from v73_count_assay import fit_log_response,invert_assay

def main():
 k=np.array([0.,.1,1.,20.,100.,10000.]);f=np.array([1000.,100.,100.,50.,5.,1.])
 effect,var=fit_log_response(k,k,f)
 np.testing.assert_array_equal(effect,np.zeros_like(k));assert np.isfinite(var).all() and (var>0).all()
 effect,_=fit_log_response(np.array([200000.,50000.]),np.array([100000.,100000.]),np.ones(2))
 np.testing.assert_allclose(effect,[1.,-1.],atol=1e-4)
 effect,_=fit_log_response(np.array([20.]),np.array([100.]),np.array([1000.]))
 assert abs(effect[0])<.5
 q=np.array([0.,1e-7,1e-5,.001,.1]);zero=np.zeros((2,len(q)))
 np.testing.assert_allclose(invert_assay(zero,q),zero,atol=1e-12)
 direction=np.array([[0.,.02,.1,.2,.3]])
 assert (invert_assay(direction,q)[0,1:]>0).all()
 print('V73 exact null under weak/high-dispersion evidence, known response, shrinkage and assay inversion passed.')

if __name__=='__main__':main()
