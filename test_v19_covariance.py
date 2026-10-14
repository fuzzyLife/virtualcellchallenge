"""Check covariance transfer, operator bounds, no-response fallback and support masks."""
import numpy as np,pandas as pd
from scipy import sparse
from v19_methods import bounded_operator,regularized_covariance,query_corrections,effects_for_target,control_stats

def main():
 np.testing.assert_allclose(bounded_operator(np.eye(16)*2,np.eye(16)),np.eye(16)*2,atol=1e-12)
 op=bounded_operator(np.eye(16)*100,np.eye(16));assert np.max(np.linalg.svd(op,compute_uv=False))<=2+1e-12
 rng=np.random.default_rng(19);x=sparse.csr_matrix(rng.poisson(rng.uniform(.5,10,50),(400,50)));depth=np.asarray(x.sum(1)).ravel();meta=pd.DataFrame({'sample':np.repeat(['batch1','batch2'],200),'gene_target':'Non-Targeting','total_counts':depth});stat=control_stats(x,meta,np.arange(400));np.testing.assert_allclose(stat['depth_residual'].reshape(2,200).mean(1),0,atol=1e-14)
 basis=np.linalg.qr(rng.normal(size=(50,16)))[0];model={'positions':np.arange(50),'scale':np.ones(50),'basis':basis,'source_covariances':np.array([np.eye(16),np.eye(16)*2]),'source_shifts':np.zeros((2,300,16)),'source_sizes':np.ones((2,300),int)*100};out,diag=query_corrections(model,x,meta,np.arange(400),np.ones(50,bool));assert all(np.count_nonzero(v)==0 for v in out.values())
 model['source_shifts'][:]=rng.normal(size=(2,300,16));out,diag=query_corrections(model,x,meta,np.arange(400),np.ones(50,bool));assert max(abs(v).max() for v in out.values())<=.250001
 measured=np.ones(50,bool);measured[:5]=False;out,_=query_corrections(model,x,meta,np.arange(400),measured);assert all(np.count_nonzero(v[:,:5])==0 for v in out.values())
 measured[:20]=False;out,diag=query_corrections(model,x,meta,np.arange(400),measured);assert diag['abstained'] and all(np.count_nonzero(v)==0 for v in out.values())
 base=rng.normal(0,.1,50).astype(np.float32);base[7]=0;effects=effects_for_target(base,out,0,np.ones(50)*100,7)
 for v in effects.values():np.testing.assert_array_equal(v,base)
 print('V19 identity map, singular bounds, nuisance centering, zero-force fallback, correction bounds, missing-gene masks, abstention and own-target preservation passed.',flush=True)
if __name__=='__main__':main()
