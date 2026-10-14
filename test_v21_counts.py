"""Check that real/shuffled joint decoders have identical architecture and bounds."""
import numpy as np
from scipy import sparse
from v21_methods import emit
from v6_methods import moments
from evaluate_atlas_transfer import load_atlas_model

def main():
 rng=np.random.default_rng(21);x=sparse.csr_matrix(rng.poisson(8,(400,50)));basis=np.linalg.qr(rng.normal(size=(50,16)))[0];q={'positions':np.arange(50),'basis':basis,'scale':np.ones(50),'derivative':np.ones(50)*.2,'query_rows':np.arange(400),'query_scores':rng.normal(size=(400,16))};q['cached_scores']=q['query_scores'].copy();program=rng.normal(0,.02,(1,400,16)).astype(np.float32);expert={'program_delta':program,'shuffled_program_delta':program.copy(),'supported':np.ones((1,50),bool),'direction':rng.normal(size=(1,16))};effect=np.zeros(50);prepared=moments(x)[:3];atlas,_=load_atlas_model();a,_=emit(x,effect,0,21,'stack_joint',prepared,atlas,expert,0,q,np.arange(400));z,_=emit(x,effect,0,21,'shuffled_joint',prepared,atlas,expert,0,q,np.arange(400));assert (a!=z).nnz==0;np.testing.assert_array_equal(np.asarray(a.sum(1)).ravel(),prepared[0]);print('V21 real/shuffled joint architecture parity, integer libraries and inherited bounds passed.',flush=True)
if __name__=='__main__':main()
