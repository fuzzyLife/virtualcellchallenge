"""Check exact cell-peer moments against explicitly enumerated Spearman pairs."""
import unittest
import numpy as np
from scipy.stats import spearmanr
from analyze_target_cell_dispersion import normalized_ranks, peer_moments, row_pair_spearman

class DispersionTests(unittest.TestCase):
    def test_exact_moments_exclude_self_and_constants(self):
        raw=np.array([[0,0,1,3,2],[2,0,0,1,4],[0,3,1,1,2],[1,1,1,1,1],[2,1,0,0,3],[0,1,0,2,4]],float)
        groups=np.array(['A','A','A','B','B','B'])
        z,valid=normalized_ranks(raw)
        results=peer_moments(z,valid,groups,batch_size=2)
        for i in range(len(raw)):
            for group,(mean,sd,count) in results.items():
                peers=[j for j in range(len(raw)) if groups[j]==group and valid[j] and i!=j]
                self.assertEqual(count[i],len(peers))
                if not valid[i] or not peers:
                    self.assertTrue(np.isnan(mean[i])); continue
                explicit=np.array([spearmanr(raw[i],raw[j]).statistic for j in peers])
                self.assertAlmostEqual(mean[i],explicit.mean(),places=13)
                self.assertAlmostEqual(sd[i],explicit.std(ddof=0),places=12)

    def test_leave_one_out_and_missing_axis(self):
        x=np.array([[0,0,1,2],[2,1,0,4],[1,0,2,3]],float)
        other=x.sum(axis=0)[None,:]-x
        actual=row_pair_spearman(x,other)
        for i in range(len(x)):
            self.assertAlmostEqual(actual[i],spearmanr(x[i],np.delete(x,i,axis=0).sum(axis=0)).statistic,places=13)
        x[0,0]=np.nan
        with self.assertRaises(ValueError):
            normalized_ranks(x)

if __name__=='__main__':
    unittest.main()
