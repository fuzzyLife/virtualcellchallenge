"""Meaningful invariant checks for the V17 generator."""
import numpy as np
from scipy import sparse
from v17_methods import emit
from v16_methods import emit as prior
from v6_methods import moments
from evaluate_atlas_transfer import load_atlas_model

def main():
 rng=np.random.default_rng(171);x=sparse.csr_matrix(rng.poisson(rng.uniform(1,20,120)[None,:],size=(400,120)).astype(np.int32));prepared=moments(x)[:3];atlas,_=load_atlas_model();effect=np.r_[np.full(30,.7),np.full(30,-.6),np.zeros(60)]
 for e,cv in [(np.zeros(120),.25),(effect,0.)]:
  a,_=emit(x,np.arange(400),e,None,17,cv,prepared,atlas);z,_=prior(x,np.arange(400),e,None,17,'residual10',prepared,atlas);assert (a!=z).nnz==0
 a,d=emit(x,np.arange(400),effect,0,17,.2,prepared,atlas);z,_=emit(x,np.arange(400),effect,0,17,.2,prepared,atlas);assert (a!=z).nnz==0
 np.testing.assert_array_equal(np.asarray(a.sum(1)).ravel(),np.asarray(x.sum(1)).ravel());assert (a.data>0).all() and np.issubdtype(a.dtype,np.integer);assert d['response_gain_variance']>0
 try:emit(x,np.arange(400),effect,None,17,.6,prepared,atlas)
 except ValueError:pass
 else:raise AssertionError('Must reject unbounded variation')
 print('V17 null identity, zero-variance V16 identity, deterministic repeat, exact libraries, integer counts, moments and bounds passed.',flush=True)
if __name__=='__main__':main()
