"""Many-small-batch depth-offset regression checks, including a known response."""
import numpy as np
from scipy import sparse
from v72_batch_offset import control_rate,response

def main():
 rng=np.random.default_rng(7272);p=np.array([.02,.05,.13,.8]);mat=[];groups=[];control={};responses=[]
 for b in range(80):
  q=p*np.array([1+.4*np.sin(b),1.,1.,1.]);q/=q.sum()
  ntc=np.array([rng.multinomial(8000,q) for _ in range(30)])
  control[str(b)]=control_rate(sparse.csr_matrix(ntc))
  # Exactly two treated cells per batch: V70's >=10 rule would discard every batch.
  for _ in range(2):mat.append(rng.multinomial(int(rng.integers(6000,20000)),q));groups.append(str(b))
 x=sparse.csr_matrix(mat,dtype=float);ids=np.arange(len(mat));groups=np.asarray(groups)
 effect,var,k,e=response(x,ids,groups,control)
 np.testing.assert_allclose(effect,np.zeros(4),atol=.035)
 assert np.isfinite(var).all() and (var>0).all() and (k>=20).all() and (e>=20).all()
 print('V72 recovered 160 treated cells across 80 small batches without depth/batch-induced response.')
 # Known multiplicative response relative to each batch's own control composition.
 changed=[]
 for b in range(80):
  q=p*np.array([1+.4*np.sin(b),1.,1.,1.]);q/=q.sum();q=q*np.array([2.,1.,1.,1.]);q/=q.sum()
  for _ in range(2):changed.append(rng.multinomial(12000,q))
 ee,_,_,_=response(sparse.csr_matrix(changed,dtype=float),ids,groups,control)
 assert .90<ee[0]<1.05 and np.max(abs(ee[1:]))<.08
 print('V72 known intervention is recovered under heterogeneous batch baselines.')

if __name__=='__main__':main()
