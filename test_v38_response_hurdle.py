"""Meaningful V38 leakage, missingness, target-mask and count-contract tests."""
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v38_response_hurdle as m
def main():
 d=m.data()
 for outer in m.OUTERS:
  ds=m.descriptor(d['profiles'][0],tuple(sorted(set([outer,'Replogle']))))
  assert all(g not in ds['excluded'] for _,g,*_ in ds['sources'])
  assert len({g for _,g,*_ in ds['sources']})==len(ds['sources'])
  ti=0;positions=np.arange(300);x,n=m.features(ds,ti,positions)
  assert x.shape==(300,56);assert not np.isinf(x).any()
  assert all(not ds['sources'][i][4][ti,ds['own'][ti]] for i in range(len(ds['sources'])))
 # Changing only an excluded study's response must not change donor features.
 saved=d['truth'].copy();first=m.descriptor(d['profiles'][0],('H1','Replogle'));xx,_=m.features(first,0,np.arange(100))
 d['truth'][np.isin(d['row_families'],['H1','K562'])]=2.9;m.bank.cache_clear()
 changed=m.descriptor(d['profiles'][0],('H1','Replogle'));yy,_=m.features(changed,0,np.arange(100));np.testing.assert_array_equal(xx,yy)
 d['truth'][:]=saved;m.bank.cache_clear()
 # Count generator retains all row libraries, nonnegative integers and own mask.
 from evaluate_atlas_transfer import load_atlas_model
 from v6_methods import moments
 atlas,_=load_atlas_model();rng=np.random.default_rng(38038);z=sparse.csr_matrix(rng.poisson(2,size=(40,300)));libs,p,bulk,_=moments(z)
 out,diag=m.emit(z,np.zeros(300),3,38038,(libs,p,bulk),atlas,False)
 np.testing.assert_array_equal(out.sum(1).A1,libs);assert out.dtype.kind in 'iu' and out.data.min(initial=0)>=0
 assert out[:,3].sum()<z[:,3].sum()*.3
 print('V38 study exclusion mutation, missingness, own-target and integer contracts PASS',flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()

