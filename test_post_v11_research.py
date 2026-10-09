"""Scientific leakage and axis invariants for the post-V11 research diagnostics."""
import numpy as np
import build_submission as build
from scipy import sparse
from audit_post_v11_transfer import fold,cosine
from benchmark_post_v11_distributions import simulate
from v7_response_methods import FAMILIES

def main():
 rng=np.random.default_rng(7);n=10;g=30;t=2
 data={'truth':rng.normal(0,.1,(n,g)).astype(np.float32),'errors':np.full((n,g),.1,np.float32),
  'valid':np.ones((n,g),bool),'row_families':np.repeat(FAMILIES,2),'row_targets':np.tile(np.arange(t),5),
  'knockdown':np.full(n,.7),'targets':np.array(['G0','G1']),'genes':np.array([f'G{i}' for i in range(g)])}
 data['valid'][0,5]=False
 base=fold(data,'K562');changed={k:v.copy() for k,v in data.items()};held=changed['row_families']=='K562'
 changed['truth'][held]+=1000;changed['errors'][held]*=100;changed['knockdown'][held]=0
 for a,b in zip(base,fold(changed,'K562')):np.testing.assert_array_equal(a,b)
 missing={k:v.copy() for k,v in data.items()};missing['truth'][~missing['valid']]=1e6
 for a,b in zip(fold(data,'HEK293T'),fold(missing,'HEK293T')):np.testing.assert_array_equal(a,b)
 zero={k:v.copy() for k,v in data.items()};zero['valid'][:]=False
 out=fold(zero,'H1');assert not out[2].any() and not out[0].any() and not out[4].any()
 x=np.array([[1.,0],[0,1.]]);np.testing.assert_allclose(cosine(x,x,np.ones_like(x,bool)),1.)
 assert np.isnan(cosine(np.zeros_like(x),x,np.ones_like(x,bool))).all()
 template=sparse.csr_matrix([[2,1,0],[0,2,2]],dtype=np.int64);depth=np.array([3,4]);p=np.array([.3,.4,.3]);var=np.array([.1,.1,.1])
 class Solver:
  def dual_moment_counts(self,profile,p,b,depths,seed,iterations):
   assert profile.shape==template.shape and np.isfinite(profile).all() and (profile>=0).all()
   np.testing.assert_allclose(profile.sum(axis=1),1.)
   return np.array([[2,1,0],[0,2,2]],dtype=np.int64)
 for method in ['gene_gamma_poisson','empirical_blend25']:
  result,_=simulate(template,depth,p,p,var,method,12,Solver());np.testing.assert_array_equal(np.asarray(result.sum(axis=1)).ravel(),depth)
 build.write_json(build.ROOT/'analysis/post_v11_research/scientific_test_results.json',{'passed':True,'script_sha256':build.sha256(build.ROOT/'test_post_v11_research.py'),'dependencies':{n:build.sha256(build.ROOT/n) for n in ['audit_post_v11_transfer.py','benchmark_post_v11_distributions.py','v7_response_methods.py']},'checks':['held-family label leakage','masked-value abstention','no-source abstention','undefined cosine handling','pilot template probability/count/depth axes'],'scope':'Scientific unit invariants; the real count-solver outputs are independently checked in 48 neutral comparison records.'})
 print('PASS: held-family response leakage, masked-value abstention, no-source abstention, undefined cosine handling, count/depth axes and valid pilot probability templates.')
if __name__=='__main__':main()
