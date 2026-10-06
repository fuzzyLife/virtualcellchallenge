"""Numerical, serialization, gradient and leakage checks for neural reconstruction."""
from argparse import Namespace
from pathlib import Path
import tempfile
import unittest
import numpy as np
from scipy import sparse
from scipy.special import erf, expit
from neural_runtime import activate_local_runtime
activate_local_runtime()
import torch
from threadpoolctl import threadpool_limits
from train_gene300_neural import (GeneNetwork, STRUCTURES, output_from_delta, save_network, load_network,
    predict_network, metadata, train_network, baseline)
from train_gene300_models import features


def independent_forward(model, inputs, weights):
    x = ((features(inputs, str(model['family'])) - model['feature_mean']) / model['feature_scale']).astype(np.float32)
    h = x
    spec = STRUCTURES[str(model['architecture'])]
    for layer, _ in enumerate(spec['hidden']):
        h = h @ weights[f'body.{layer * 3}.weight'].T + weights[f'body.{layer * 3}.bias']
        if spec['activation'] == 'relu':
            h = np.maximum(h, 0)
        elif spec['activation'] == 'silu':
            h = h * expit(h)
        else:
            h = h * .5 * (1 + erf(h / np.sqrt(2)))
    delta = h @ weights['head.weight'].T + weights['head.bias']
    if str(model['architecture']) == 'residual':
        delta += (x @ weights['branch.0.weight'].T) @ weights['branch.1.weight'].T
    base = baseline(model, inputs)
    if str(model['objective']) == 'raw':
        return base * np.exp(np.clip(delta, -3, 3))
    return np.clip(base + delta, 0, np.log1p(1e6))


class NeuralTests(unittest.TestCase):
    def setUp(self):
        torch.manual_seed(17)
        self.rng = np.random.default_rng(17)

    def test_all_structures_dimensions_initial_baseline_and_gradients(self):
        for arch in STRUCTURES:
            network = GeneNetwork(arch, 302, 9)
            x = torch.randn(12, 302)
            delta = network(x)
            self.assertEqual(tuple(delta.shape), (12, 9))
            torch.testing.assert_close(delta, torch.zeros_like(delta))
            (delta - 1).square().mean().backward()
            self.assertGreater(network.head.bias.grad.abs().sum().item(), 0)
            self.assertTrue(all(torch.isfinite(p.grad).all() for p in network.parameters() if p.grad is not None))

    def test_raw_head_is_positive_and_zero_baseline_stays_zero(self):
        base = torch.tensor([[0., 1., 100.], [2., 4., 10.]])
        delta = torch.tensor([[-100., 0., 100.], [-1., 1., 0.]])
        out = output_from_delta(delta, base, 'raw')
        self.assertTrue((out >= 0).all() and torch.isfinite(out).all())
        self.assertEqual(out[0, 0].item(), 0)
        self.assertAlmostEqual(out[0, 2].item(), 100 * np.exp(3), places=3)

    def test_both_objectives_metadata_ignore_heldout_outputs(self):
        inputs = self.rng.poisson(3, (12, 300)).astype(float)
        raw = self.rng.poisson(4, (12, 7))
        labels = np.log1p(raw * 100).astype(np.float32)
        fit = np.arange(12) < 8
        changed = raw.copy(); changed[~fit] = 99999
        ychanged = labels.copy(); ychanged[~fit] = 999
        for objective in ('normalized', 'raw'):
            a = metadata(inputs, sparse.csr_matrix(raw), labels, inputs.sum(axis=1), fit, 'A', 'deep', objective, 'test', 5)
            b = metadata(inputs, sparse.csr_matrix(changed), ychanged, inputs.sum(axis=1), fit, 'A', 'deep', objective, 'test', 5)
            for key in a:
                np.testing.assert_array_equal(a[key], b[key])

    def test_checkpoint_replay_and_independent_numpy_forward_all_structures(self):
        inputs = self.rng.poisson(3, (13, 300)).astype(float)
        for objective in ('normalized', 'raw'):
            for arch in STRUCTURES:
                family = STRUCTURES[arch]['family']
                f = features(inputs, family)
                model = {'schema': np.array(1), 'architecture': np.array(arch), 'objective': np.array(objective),
                    'family': np.array(family), 'feature_mean': f.mean(axis=0),
                    'feature_scale': np.maximum(f.std(axis=0), 1e-6), 'base_vector': np.ones(9),
                    'input_genes': np.array([f'g{i}' for i in range(300)]), 'output_genes': np.array([f'y{i}' for i in range(9)])}
                network = GeneNetwork(arch, f.shape[1], 9).eval()
                with torch.no_grad():
                    network.head.weight.normal_(0, .05)
                    if network.branch is not None:
                        network.branch[-1].weight.normal_(0, .01)
                actual = predict_network(model, network, inputs)
                expected = independent_forward(model, inputs, {k: v.detach().numpy() for k, v in network.state_dict().items()})
                np.testing.assert_allclose(actual, expected, rtol=1e-5, atol=1e-4)
                with tempfile.TemporaryDirectory() as folder:
                    path = Path(folder) / 'model.npz'
                    save_network(path, model, network)
                    loaded, net = load_network(path)
                    np.testing.assert_allclose(predict_network(loaded, net, inputs), actual, rtol=1e-6, atol=1e-6)

    def test_early_stopping_and_test_label_isolation(self):
        inputs = self.rng.poisson(3, (16, 300)).astype(float)
        y = np.log1p(self.rng.poisson(3, (16, 5))).astype(np.float32)
        fit, val = np.arange(16) < 10, (np.arange(16) >= 10) & (np.arange(16) < 13)
        template = metadata(inputs, sparse.csr_matrix(np.expm1(y).astype(int)), y, inputs.sum(axis=1), fit,
                            'A', 'bottleneck', 'normalized', 'test', 5)
        template['input_genes'] = np.array([f'g{i}' for i in range(300)])
        template['output_genes'] = np.array([f'y{i}' for i in range(5)])
        args = Namespace(device='cpu', learning_rate=.001, epochs=3, patience=2, batch_size=8)
        changed = y.copy(); changed[~fit & ~val] = 999
        models, curves, metas = [], [], []
        for labels in (y, changed):
            model = dict(template)
            net, curve = train_network(model, inputs, labels, fit, val, np.ones(5, bool), args, 'A', 'bottleneck', 'normalized')
            models.append(net); curves.append(curve); metas.append(model)
            self.assertAlmostEqual(float(model['validation_mse']), min(r['validation_mse'] for r in curve))
        self.assertEqual(curves[0], curves[1])
        for key in models[0].state_dict():
            torch.testing.assert_close(models[0].state_dict()[key], models[1].state_dict()[key], rtol=0, atol=0)


if __name__ == '__main__':
    torch.set_num_threads(2)
    with threadpool_limits(limits=2):
        unittest.main(verbosity=2)
