"""Focused numerical and leakage-boundary checks for the new reconstruction models."""
import tempfile
from pathlib import Path
import unittest
import numpy as np
from scipy.linalg import lstsq
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits
from train_gene300_models import (features, standardize, split_guides, ridge_fit, ridge_basis,
    validation_sse, fit_model, predict, load_model, spearman)
from predict_gene300_models import read_inputs


class GeneModelTests(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(14)

    def test_guide_splits(self):
        guides = np.repeat([f'g{i:02}' for i in range(46)], 400)
        roles = split_guides(guides)
        np.testing.assert_array_equal(roles, split_guides(guides))
        self.assertEqual({r: int((roles == r).sum()) for r in set(roles)},
                         {'train': 12800, 'validation': 2800, 'test': 2800})
        self.assertTrue(all(len(set(roles[guides == g])) == 1 for g in np.unique(guides)))

    def test_zero_ties_and_rank_order(self):
        x = np.tile(np.array([0., 0., 1., 3., 3., 7.]), 50)[None, :]
        a, b = features(x, 'rank'), features(x ** 2, 'rank')
        np.testing.assert_allclose(a[:, :300], b[:, :300])
        self.assertEqual(a[0, 0], a[0, 1])
        np.testing.assert_array_equal(features(np.zeros((1, 300)), 'rank'), np.zeros((1, 302)))
        for bad in (np.full((1, 300), np.nan), np.full((1, 300), -1), np.zeros((1, 301))):
            with self.assertRaises(ValueError):
                features(bad, 'hybrid')

    def test_ridge_against_augmented_least_squares(self):
        x, _, _ = standardize(self.rng.normal(size=(80, 7)))
        y = self.rng.normal(size=(80, 4))
        alpha = .3
        coef, intercept = ridge_fit(x, y, alpha)
        expected = lstsq(np.vstack((x, np.sqrt(len(x) * alpha) * np.eye(7))),
                         np.vstack((y - y.mean(axis=0), np.zeros((7, 4)))))[0]
        np.testing.assert_allclose(coef, expected, atol=1e-12)
        np.testing.assert_allclose(intercept, y.mean(axis=0))

    def test_sufficient_statistics_match_explicit_loss(self):
        x, _, _ = standardize(self.rng.normal(size=(80, 7)))
        xv, y, yv = self.rng.normal(size=(20, 7)), self.rng.normal(size=(80, 4)), self.rng.normal(size=(20, 4))
        values, vectors, z = ridge_basis(x)
        zv = xv @ vectors
        mean = y.mean(axis=0)
        b = (z.T @ (y - mean) / len(x)) / (values[:, None] + .1)
        got = validation_sse(b, zv.T @ (yv - mean), zv.T @ zv, ((yv - mean) ** 2).sum(axis=0))
        np.testing.assert_allclose(got, ((zv @ b + mean - yv) ** 2).sum(axis=0), atol=1e-12)

    def test_training_scope_serialization_and_count_conservation(self):
        inputs = self.rng.poisson(2, size=(70, 300)).astype(float)
        labels = np.log1p(self.rng.poisson(3, size=(70, 12)) * 100.)
        library = inputs.sum(axis=1) + self.rng.integers(1000, 2000, size=70)
        mask = np.arange(70) < 50
        selected = {'family': 'magnitude', 'alpha': 1.}
        names = np.array([f'in{i}' for i in range(300)], dtype=object)
        model = fit_model(inputs, labels, library, mask, selected, names, [f'out{i}' for i in range(12)], 'A', 'test', 1.)
        # Changing every held-out label/depth leaves fitted weights and statistics identical.
        labels2, lib2 = labels.copy(), library.copy()
        labels2[~mask] = 1000; lib2[~mask] = 999999
        other = fit_model(inputs, labels2, lib2, mask, selected, names, [f'out{i}' for i in range(12)], 'A', 'test', 1.)
        for key in model:
            np.testing.assert_array_equal(model[key], other[key])
        counts, latent, depth = predict(model, inputs[~mask])
        self.assertTrue(np.isfinite(counts).all() and (counts >= 0).all())
        np.testing.assert_allclose(counts.sum(axis=1, dtype=np.float64) + inputs[~mask].sum(axis=1), depth, rtol=1e-7)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'model.npz'
            np.savez_compressed(path, **model)
            reloaded = load_model(path)
            np.testing.assert_allclose(predict(reloaded, inputs[~mask])[0], counts, rtol=1e-7)
        batched = np.vstack([predict(model, row[None, :])[0] for row in inputs[~mask]])
        np.testing.assert_allclose(batched, counts, rtol=1e-6, atol=1e-6)

    def test_correlations_keep_zeros_and_undefined_constants(self):
        x = np.array([0, 0, 0, 1, 3, 4])[:, None]
        y = np.array([0, 2, 0, 3, 1, 5])[:, None]
        self.assertAlmostEqual(spearman(x, y)[0], spearmanr(x[:, 0], y[:, 0]).statistic)
        self.assertTrue(np.isnan(spearman(x, np.ones_like(x))[0]))

    def test_named_input_reordering_and_validation(self):
        names = np.array([f'g{i}' for i in range(300)])
        import pandas as pd
        frame = pd.DataFrame(np.arange(600).reshape(2, 300), columns=names)
        frame['context'] = ['A', 'B']
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'input.csv'
            frame[list(frame.columns)[::-1]].to_csv(path, index=False)
            inputs, obs = read_inputs(path, {'input_genes': names}, 'A')
            np.testing.assert_array_equal(inputs[0], np.arange(300))
            self.assertEqual(len(obs), 1)
            frame.drop(columns=['g1']).to_csv(path, index=False)
            with self.assertRaises(ValueError):
                read_inputs(path, {'input_genes': names}, 'A')


if __name__ == '__main__':
    with threadpool_limits(limits=2):
        unittest.main(verbosity=2)
