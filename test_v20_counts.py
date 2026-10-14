"""Scientific invariants for Stack baseline calibration and conditional count emission."""
import numpy as np
from scipy import sparse
from v20_methods import summarize_pair,target_effects,emit,conditional_rows
from v16_methods import emit as prior
from v6_methods import moments
from evaluate_atlas_transfer import load_atlas_model

def main():
 rng=np.random.default_rng(20);genes=np.asarray(['g'+str(i) for i in range(50)]);basis=np.linalg.qr(rng.normal(size=(50,16)))[0];mu=rng.uniform(.1,10,(400,50));q={'positions':np.arange(50),'scale':np.ones(50),'basis':basis};summary=summarize_pair(mu,mu,mu.mean(0),q,np.ones(50,bool),np.ones(50,bool),list(genes),genes,np.ones_like(mu),np.ones_like(mu));assert np.count_nonzero(summary['lfc'])==np.count_nonzero(summary['shrink'])==np.count_nonzero(summary['program_delta'])==np.count_nonzero(summary['direction'])==0
 x=sparse.csr_matrix(rng.poisson(rng.uniform(1,20,50),(400,50)));prepared=moments(x)[:3];atlas,_=load_atlas_model();base=rng.normal(0,.1,50).astype(np.float32);base[3]=0;expert={k:np.zeros((300,50),np.float32) for k in ['lfc','shrink','raw','shuffled_lfc','shuffled_shrink']};expert.update(supported=np.zeros((300,50),bool),program_delta=np.zeros((300,400,16),np.float32),direction=np.zeros((300,16),np.float32));qi={'query_rows':np.arange(400),'query_scores':rng.normal(size=(400,16)),'positions':np.arange(50),'scale':np.ones(50),'basis':basis,'derivative':np.ones(50)*.2};qi['cached_scores']=qi['query_scores'].copy();score,nearest=conditional_rows(qi,np.arange(400));np.testing.assert_array_equal(nearest,np.arange(400))
 for kind in ['stack_joint','reweight']:
  actual,_=emit(x,base,3,20,kind,prepared,atlas,expert,0,qi,np.arange(400));expected,_=prior(x,np.arange(400),base,3,20,'residual10',prepared,atlas);assert (actual!=expected).nnz==0
 effects=target_effects(base,expert,0,np.ones(50)*100,3)
 for kind in ['v16','stack_mean','stack_joint','raw_mean','shuffled','reweight','amplitude']:np.testing.assert_array_equal(effects[kind],base)
 expert['supported'][0]=True;expert['program_delta'][0]=rng.normal(0,.02,(400,16));expert['direction'][0]=rng.normal(size=16)*.2
 for kind in ['stack_joint','reweight']:
  actual,diag=emit(x,base,3,20,kind,prepared,atlas,expert,0,qi,np.arange(400));repeat,_=emit(x,base,3,20,kind,prepared,atlas,expert,0,qi,np.arange(400));assert (actual!=repeat).nnz==0;np.testing.assert_array_equal(np.asarray(actual.sum(1)).ravel(),prepared[0]);assert (actual.data>0).all() and np.issubdtype(actual.dtype,np.integer)
 expert['lfc'][0]=10;expert['shrink'][0]=1;effect=target_effects(base,expert,0,np.ones(50)*100,3)['stack_mean'];assert abs(effect-base).max()<=.125001 and effect[3]==0
 print('V20 paired-control bias cancellation, missing-expert fallback, nearest-state identity, exact zero-program V16 identity, integer libraries, seeded conditional replay, own-target exclusion and bounded updates passed.',flush=True)
if __name__=='__main__':main()
