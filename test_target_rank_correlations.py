"""Independent oracles for pairwise-valid raw rank correlations."""
import unittest
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from analyze_target_rank_correlations import pairwise_spearman, pairwise_rank_rows

class TargetRankTests(unittest.TestCase):
    def test_pairwise_missing_zeros_and_ties(self):
        cells = np.array([[0, 0, 2, 4, 1, 9], [3, np.nan, 0, 0, 5, 4], [np.nan, 1, 3, 2, 0, 2]], float)
        reference = np.array([0, 3, np.nan, 4, 1, 1], float)
        actual, counts, statuses = pairwise_rank_rows(cells, reference)
        for i, row in enumerate(cells):
            mask = np.isfinite(row) & np.isfinite(reference)
            expected = spearmanr(row[mask], reference[mask]).statistic
            self.assertAlmostEqual(actual[i], expected, places=13)
            self.assertAlmostEqual(actual[i], pd.Series(row).corr(pd.Series(reference), method='spearman'), places=13)
            self.assertEqual(counts[i], mask.sum())
            self.assertEqual(statuses[i], 'ok')
        self.assertEqual(counts[0], 5)  # Includes both measured zero observations.

    def test_rerank_after_pairwise_removal(self):
        left = np.array([1, 50, 2, 3, 4], float)
        right = np.array([9, np.nan, 1, np.inf, 5], float)
        rho, count, status = pairwise_spearman(left, right)
        self.assertEqual(count, 3)
        self.assertEqual(status, 'ok')
        self.assertAlmostEqual(rho, spearmanr([1, 2, 4], [9, 1, 5]).statistic, places=13)

    def test_undefined_and_scale_invariance(self):
        self.assertEqual(pairwise_spearman([0, 0, 0], [1, 2, 3])[2], 'constant_rank_vector')
        self.assertTrue(np.isnan(pairwise_spearman([0, 0, 0], [1, 2, 3])[0]))
        self.assertEqual(pairwise_spearman([np.nan, 2], [1, 3])[2], 'insufficient_pairs')
        self.assertEqual(pairwise_spearman([np.inf, 2], [1, 3])[1], 1)
        x = np.array([[0, 0, 2, 7], [1, 4, 1, 8]], float)
        reference = np.array([0, 1, 3, 4], float)
        np.testing.assert_allclose(pairwise_rank_rows(x, reference)[0], pairwise_rank_rows(x*123, reference*8)[0])
        rho, count, status = pairwise_rank_rows(x, np.full(4, np.nan))
        self.assertTrue(np.isnan(rho).all())
        self.assertTrue((count == 0).all())
        self.assertTrue((status == 'insufficient_pairs').all())

if __name__ == '__main__':
    unittest.main()
