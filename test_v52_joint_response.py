"""V52 boundary and independent rank-label checks before fitting."""
import numpy as np,torch
from scipy import stats,sparse
import v52_joint_response as m

def main():
 rng=np.random.default_rng(52)
 for i in range(12):
  x=rng.integers(0,8,(17,13));c=rng.integers(0,8,(23,13));theta=m.rank_effect(x,c)
  ref=stats.mannwhitneyu(x,c,axis=0,method='asymptotic',use_continuity=False).statistic/(17*23)-.5
  np.testing.assert_allclose(theta,ref,rtol=0,atol=1e-14)
 np.testing.assert_array_equal(m.rank_effect(np.zeros((3,2)),np.zeros((7,2))),[0,0])
 # The own gene is absent from the response labels; unchanged counts rank at zero.
 x=sparse.csr_matrix(np.tile([5,2,1],(32,1)))
 mean,rank,valid,q,n=m.summarize_views(x,[[(list(range(8)),list(range(8,16)))],[(list(range(16,24)),list(range(24,32)))]],np.ones(3,bool),1)
 np.testing.assert_allclose(mean,0);np.testing.assert_allclose(rank,0);assert not valid[1]
 assert m.GROUP['HCT116']==m.GROUP['HEK293T'] and all(m.GROUP[s]=='CD4' for s in m.SOURCES if s.startswith('CD4'))
 assert not any('D3' in s for s in m.SOURCES)
 # Model bounds and the mean/rank single-head ablations are enforced in forward.
 torch.manual_seed(1);net=m.Network(10);net.eval();xx=torch.randn(7,10);prior=torch.tensor(np.column_stack([rng.uniform(-1,1,7),rng.uniform(-.1,.1,7)]),dtype=torch.float32)
 with torch.no_grad():
  p=net(xx,prior,'rank_only');np.testing.assert_allclose(p[:,0],prior[:,0]);p=net(xx,prior,'mean_only');np.testing.assert_allclose(p[:,1],prior[:,1])
 # Missingness survives imputation rather than being encoded as a biological zero.
 a=np.array([[1.,np.nan],[2.,3.],[4.,5.]],np.float32);med,scale=m.normalizer(a);v=m.transform(a,med,scale);assert v.shape==(3,4) and v[0,3]==1
 print('V52 tests passed: 156 independent U comparisons, ties, null, own-gene mask, study grouping, ablations and missingness',flush=True)
if __name__=='__main__':main()
