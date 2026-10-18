"""Numerical contracts for V59, including control-only state calibration."""
import numpy as np
from scipy import sparse
from v59_state_composition import estimate,reconstruct,standardize,logratio

def main():
 cm=np.array([[1.,10.],[10.,1.],[2.,2.],[2.,2.]])
 n=np.array([40.,40.,10.,10.]);l,w,v=estimate(n,n,cm,cm)
 np.testing.assert_array_equal(l,np.zeros(4));np.testing.assert_array_equal(w,np.zeros_like(cm))
 np.testing.assert_allclose(reconstruct(n,cm,l,w),(n/n.sum())@cm)
 out=reconstruct(n,cm,np.array([1.,-1.,0.,0.]));assert out[0]<((n/n.sum())@cm)[0] and out[1]>((n/n.sum())@cm)[1]
 nn=np.array([0.,40.,10.,10.]);_,_,valid=estimate(nn,n,cm,cm);assert not valid[0].any()
 m={'lane':np.array(['a']*6),'label':np.array(['Non-Targeting']*4+['X']*2),'rows':np.arange(6)}
 x=sparse.csr_matrix(np.arange(1,19).reshape(6,3));z,t=standardize(x,m,np.arange(3));xx=x.toarray();xx[4:]*=100;zz,tt=standardize(sparse.csr_matrix(xx),m,np.arange(3));assert t==tt;np.testing.assert_array_equal(z[:4],zz[:4])
 assert np.isfinite(logratio(np.zeros(3),np.zeros(3))).all()
 print('V59 numerical/control-only contracts passed',flush=True)
if __name__=='__main__':main()

# Count editing must conserve each original integer library and preserve zero intervention.
def count_contract():
 from evaluate_v59_counts import state_counts
 rng=np.random.default_rng(59);x=sparse.csr_matrix(rng.poisson(8,size=(500,3)));pool=np.arange(500);state=np.arange(500)%4
 count,rows=state_counts(x,pool,state,np.zeros(4),np.zeros((4,3)),np.zeros(3),None,59)
 assert (count!=x[rows]).nnz==0
 count,rows=state_counts(x,pool,state,np.array([.2,-.2,0,0]),np.tile([.2,-.2,0],(4,1)),np.array([.1,0,-.1]),0,59)
 assert count.shape==(400,3) and len(set(rows))==400 and (count.data>=0).all()
 np.testing.assert_array_equal(count.sum(1),x[rows].sum(1))
 print('V59 integer-library and zero-intervention count contracts passed',flush=True)
if __name__=='__main__':count_contract()
