"""Depth/exposure and uncertainty checks for the V70 count teacher."""
import numpy as np
from scipy import sparse
from v70_count_teacher import rate

def main():
 rng=np.random.default_rng(7070)
 depth=np.repeat([2000.,5000.,10000.,20000.],200)
 p=np.array([.01,.03,.2,.76]);x=rng.multinomial(1,p,size=1)
 control=sparse.csr_matrix(np.array([rng.multinomial(int(d),p) for d in depth]))
 treated=sparse.csr_matrix(np.array([rng.multinomial(int(d*2),p) for d in depth]))
 a,av,_,_=rate(control);b,bv,_,_=rate(treated)
 np.testing.assert_allclose((b-a)/np.log(2),np.zeros(4),atol=.025)
 assert (av>0).all() and (bv>0).all() and np.isfinite(a).all()
 factor=np.array([2.,1.,1.,1.]);pt=p*factor;pt/=pt.sum()
 intervention=sparse.csr_matrix(np.array([rng.multinomial(int(d),pt) for d in depth]))
 c,cv,_,_=rate(intervention)
 np.testing.assert_allclose((c-a)/np.log(2),np.log2(pt/p),atol=.035)
 # Repeating counts raises estimated precision: report this as technical, not biological replication.
 _,small,_,_=rate(control[:40]);assert np.median(av)<np.median(small)
 zero=sparse.hstack([control,sparse.csr_matrix((len(depth),1))],format='csr')
 zz,zv,_,_=rate(zero);assert np.isfinite(zz).all() and np.isfinite(zv).all()
 print('V70 depth-offset, known-effect, uncertainty and zero-count checks passed.')

if __name__=='__main__':main()
