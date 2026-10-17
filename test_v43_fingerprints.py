"""Boundary tests for experimental-feature transfer, separate from score selection."""
import numpy as np
import v43_lincs_fingerprints as m
def main():
 keys,e=m.embeddings();d=m.data()
 assert e.shape==(len(keys),64) and len(set(keys))==len(keys) and np.isfinite(e).all()
 np.testing.assert_allclose(np.linalg.norm(e,axis=1),1,atol=1e-6)
 with np.load(m.A/'experimental_fingerprints.npz',allow_pickle=False) as z:
  assert not np.isin(z['landmarks'],d['targets']).any();assert (z['context_support']>=3).all()
 m.install()
 # Numerically exact own-coordinate exclusion in the reused solver.
 rng=np.random.default_rng(43);x=rng.normal(size=(20,4));g=rng.normal(size=(12,3));y=rng.normal(size=(20,12));own=np.arange(20)%12;w=np.ones(20)
 a=m.core.masked_bilinear(x,g,y,own,w,2.);y[np.arange(20),own]+=100000
 b=m.core.masked_bilinear(x,g,y,own,w,2.);np.testing.assert_allclose(a,b,atol=1e-9)
 m.write(m.A/'contract_tests.json',{'passed':True,'checks':['unique finite 64D experimental embeddings','unit L2 norms','all300 target outputs absent from fingerprints','>=3 source contexts','own labels cannot change fitted bilinear weights'],'script_sha256':m.b.sha256(m.ROOT/'test_v43_fingerprints.py')})
 print('V43 five contract checks passed',flush=True)
if __name__=='__main__':main()

