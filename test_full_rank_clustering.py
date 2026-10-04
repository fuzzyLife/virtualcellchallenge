import unittest
import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score
from cluster_full_rank_matrix import rank_unit_rows, spherical_refine


class FullRankClusteringTests(unittest.TestCase):
    def test_tied_zero_rank_geometry(self):
        raw = np.array([[0, 0, 4, 1, 0, 2], [0, 2, 1, 1, 0, 8], [8, 8, 4, 1, 0, 2]])
        z = rank_unit_rows(raw)
        for i in range(len(raw)):
            for j in range(len(raw)):
                self.assertAlmostEqual(float(z[i] @ z[j]), spearmanr(raw[i], raw[j]).statistic, places=6)
        np.testing.assert_allclose(np.sum(z, axis=1), 0, atol=1e-7)
        np.testing.assert_allclose(np.sum(z*z, axis=1), 1, atol=1e-7)

    def test_missing_and_constant_vectors_are_rejected(self):
        for raw in [[[0, np.nan, 1]], [[0, np.inf, 1]], [[0, 0, 0]]]:
            with self.assertRaises(ValueError):
                rank_unit_rows(raw)

    def test_blind_spherical_fit_recovers_distinct_rank_patterns(self):
        raw = np.array([[0, 0, 1, 3, 8, 20], [0, 0, 2, 4, 9, 21],
                        [20, 8, 3, 1, 0, 0], [21, 9, 4, 2, 0, 0],
                        [0, 20, 1, 8, 3, 0], [0, 21, 2, 9, 4, 0]])
        z = rank_unit_rows(raw)
        labels, centers, history, converged = spherical_refine(z, z[[0, 2, 4]])
        self.assertTrue(converged)
        self.assertEqual(adjusted_rand_score([0, 0, 1, 1, 2, 2], labels), 1)
        np.testing.assert_array_equal(labels, (z @ centers.T).argmax(axis=1))
        self.assertTrue(np.all(np.diff([x['sum_centroid_rho'] for x in history]) >= -1e-5))
        np.testing.assert_allclose(np.linalg.norm(centers, axis=1), 1, atol=1e-6)
        for k in range(3):
            mean = z[labels == k].mean(axis=0)
            np.testing.assert_allclose(centers[k], mean / np.linalg.norm(mean), atol=1e-6)


if __name__ == '__main__':
    unittest.main()
