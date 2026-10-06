"""Independently verify neural models, frozen selection, metrics and output files."""
import argparse
import json
from pathlib import Path
from neural_runtime import activate_local_runtime
activate_local_runtime()
import anndata as ad
import h5py
import numpy as np
import pandas as pd
import torch
from threadpoolctl import threadpool_limits
from build_submission import ROOT, read_spec, sha256, write_json
from train_gene300_models import features, split_guides
from train_gene300_neural import OUT, RIDGE, load_network, predict_network
from validate_gene300_models import strings
from test_gene300_neural import independent_forward


def audit(out):
    manifest, genes, targets = read_spec(ROOT)
    provenance = json.loads((out / 'run_manifest.json').read_text(encoding='utf-8'))
    assert set(provenance['contexts']) == set('ABC')
    assert sha256(ROOT / 'train_gene300_neural.py') == provenance['code_sha256']
    for name, record in provenance['outputs'].items():
        assert (out / name).stat().st_size == record['bytes']
        assert sha256(out / name) == record['sha256'], name
    candidates = pd.read_csv(out / 'candidate_scores.csv')
    scores = pd.read_csv(out / 'test_architecture_scores.csv')
    per_gene = pd.read_csv(out / 'per_gene_test_metrics.csv')
    per_cell = pd.read_csv(out / 'per_cell_test_metrics.csv')
    assert len(candidates) == 24 and len(scores) == 24
    assert len(per_gene) == 109398 and len(per_cell) == 16800
    ipos = [np.flatnonzero(genes == gene)[0] for gene in targets]
    opos = [i for i, gene in enumerate(genes) if gene not in set(targets)]
    report = {'passed': True, 'source_and_output_hashes_match': True,
        'architecture_candidates': 24, 'contexts': {},
        'independent_math': 'NumPy forward pass with explicit ReLU/SiLU/GELU and low-rank linear branch',
        'prediction_coverage': 'all held-out and full raw values scanned; 32 rows per model independently replayed'}
    for context in 'ABC':
        print(context, 'checking neural models and every stored prediction value', flush=True)
        assert sha256(ROOT / f'context_{context}.h5ad') == provenance['contexts'][context]['source_sha256']
        source = ad.read_h5ad(ROOT / f'context_{context}.h5ad')
        raw = source.X.astype(np.int32).tocsr(); raw.sum_duplicates(); raw.eliminate_zeros()
        panel = raw[:, ipos].toarray().astype(np.float64)
        libraries = np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel()
        roles = split_guides(source.obs.ntc_id.astype(str).to_numpy())
        train, test = roles == 'train', roles == 'test'
        test_rows = np.flatnonzero(test)
        selected = json.loads((out / f'selected_{context}.json').read_text(encoding='utf-8'))
        for objective in ('normalized', 'raw'):
            selection = selected[objective]
            valid = candidates[(candidates.context == context) & (candidates.objective == objective)]
            winner = valid.sort_values('validation_mse').iloc[0]
            assert selection['architecture'] == winner.architecture
            assert selection['best_epoch'] == winner.best_epoch
            for prefix, fit, scope in [('evaluation', train, 'train_only_evaluation'),
                                       ('model', np.ones(len(panel), dtype=bool), 'all_controls_production')]:
                model, network = load_network(out / f'{prefix}_{objective}_{context}.npz')
                assert str(model['context']) == context and str(model['objective']) == objective
                assert str(model['architecture']) == selection['architecture'] and str(model['scope']) == scope
                assert int(model['best_epoch']) == selection['best_epoch']
                assert int(model['seed']) == selection['seed']
                np.testing.assert_array_equal(model['fit_row_indices'], np.flatnonzero(fit))
                np.testing.assert_array_equal(model['input_genes'], targets)
                np.testing.assert_array_equal(model['output_genes'], genes[opos])
                f = features(panel[fit], str(model['family']))
                scale = f.std(axis=0); scale[scale < 1e-12] = 1
                np.testing.assert_allclose(model['feature_mean'], f.mean(axis=0), atol=1e-12)
                np.testing.assert_allclose(model['feature_scale'], scale, atol=1e-12)
                if objective == 'raw':
                    expected = np.asarray(raw[fit][:, opos].sum(axis=0)).ravel() / panel[fit].sum()
                    np.testing.assert_allclose(model['base_vector'], expected, rtol=1e-12)
                elif 'ridge_coefficient' in model:
                    ridge_file = 'evaluation_model' if prefix == 'evaluation' else 'model'
                    with np.load(RIDGE / f'{ridge_file}_{context}.npz', allow_pickle=False) as ridge:
                        for key in ['family', 'feature_mean', 'feature_scale', 'coefficient', 'intercept']:
                            np.testing.assert_array_equal(model['ridge_' + key], ridge[key])
                elif objective == 'normalized':
                    sampled = np.linspace(0, 18232, 32, dtype=int)
                    counts = raw[fit][:, np.asarray(opos)[sampled]].toarray().astype(np.float64)
                    values = np.log1p(counts * (1e6 / libraries[fit])[:, None]).astype(np.float32)
                    np.testing.assert_allclose(model['base_vector'][sampled], values.mean(axis=0, dtype=np.float64), atol=1e-12)
                rows = np.linspace(0, 18399, 32, dtype=int)
                weights = {k: v.detach().numpy() for k, v in network.state_dict().items()}
                independent = independent_forward(model, panel[rows], weights)
                np.testing.assert_allclose(predict_network(model, network, panel[rows]), independent, rtol=3e-5, atol=2e-5)
            model, network = load_network(out / f'evaluation_{objective}_{context}.npz')
            gene_metrics = per_gene[(per_gene.context == context) & (per_gene.objective == objective)]
            cell_metrics = per_cell[(per_cell.context == context) & (per_cell.objective == objective)]
            np.testing.assert_array_equal(gene_metrics.gene_name, genes[opos])
            np.testing.assert_array_equal(cell_metrics.cell_id, source.obs_names[test])
            squared = np.zeros(18233); base_squared = squared.copy()
            cell_squared = np.zeros(2800)
            expressed = []
            with h5py.File(out / f'heldout_{objective}_{context}.h5ad', 'r') as handle:
                assert handle['X'].shape == (2800, 18233)
                np.testing.assert_array_equal(strings(handle['var/gene_name']), genes[opos])
                np.testing.assert_array_equal(strings(handle['obs/_index']), source.obs_names[test])
                np.testing.assert_array_equal(handle['obsm/measured_300_raw_counts'][:], panel[test])
                predictions = handle['X'][:].astype(np.float64)
                assert np.isfinite(predictions).all() and (predictions >= 0).all()
                for start in range(0, 18233, 256):
                    stop = min(start + 256, 18233)
                    train_raw = raw[train][:, opos[start:stop]].toarray().astype(np.float64)
                    actual = raw[test][:, opos[start:stop]].toarray().astype(np.float64)
                    train_cpm = train_raw * (1e6 / libraries[train])[:, None]
                    expressed.extend(train_cpm.mean(axis=0) >= 5)
                    if objective == 'normalized':
                        actual = np.log1p(actual * (1e6 / libraries[test])[:, None]).astype(np.float32).astype(np.float64)
                        reference = np.log1p(train_cpm).astype(np.float32).mean(axis=0, dtype=np.float64)[None, :]
                    else:
                        reference = panel[test].sum(axis=1)[:, None] * train_raw.sum(axis=0)[None, :] / panel[train].sum()
                    error = (predictions[:, start:stop] - actual) ** 2
                    squared[start:stop] = error.sum(axis=0)
                    base_squared[start:stop] = ((reference - actual) ** 2).sum(axis=0)
                    cell_squared += error.sum(axis=1)
                np.testing.assert_allclose(gene_metrics.test_mse, squared / 2800, rtol=1e-9, atol=1e-10)
                np.testing.assert_allclose(gene_metrics.baseline_test_mse, base_squared / 2800, rtol=1e-9, atol=1e-10)
                np.testing.assert_array_equal(gene_metrics.expressed_train_5cpm, expressed)
                np.testing.assert_allclose(cell_metrics.test_mse, cell_squared / 18233, rtol=1e-9, atol=1e-10)
                mask = np.asarray(expressed) if objective == 'normalized' else np.ones(18233, bool)
                score = scores[(scores.context == context) & (scores.objective == objective) & scores.selected].iloc[0]
                np.testing.assert_allclose(score.test_mse, squared[mask].mean() / 2800, rtol=1e-9)
                np.testing.assert_allclose(score.test_skill, 1 - squared[mask].sum() / base_squared[mask].sum(), atol=1e-12)
                rows = np.linspace(0, 2799, 32, dtype=int)
                independent = independent_forward(model, panel[test_rows[rows]], {k: v.detach().numpy() for k, v in network.state_dict().items()})
                np.testing.assert_allclose(predictions[rows], independent, rtol=3e-5, atol=2e-5)
                if objective == 'raw':
                    np.testing.assert_allclose(handle['obs/predicted_library_counts'][:], predictions.sum(axis=1) + panel[test].sum(axis=1), rtol=1e-6)
            del predictions
        model, network = load_network(out / f'model_raw_{context}.npz')
        with h5py.File(out / f'reconstructed_raw_{context}.h5ad', 'r') as handle:
            assert handle['X'].shape == (18400, 18233)
            np.testing.assert_array_equal(strings(handle['var/gene_name']), genes[opos])
            np.testing.assert_array_equal(strings(handle['obs/_index']), source.obs_names)
            np.testing.assert_array_equal(handle['obsm/measured_300_raw_counts'][:], panel)
            for start in range(0, 18400, 256):
                values = handle['X'][start:start + 256]
                assert np.isfinite(values).all() and (values >= 0).all()
                np.testing.assert_allclose(values.sum(axis=1, dtype=np.float64) + panel[start:start + 256].sum(axis=1),
                                           handle['obs/predicted_library_counts'][start:start + 256], rtol=1e-6)
            rows = np.linspace(0, 18399, 32, dtype=int)
            independent = independent_forward(model, panel[rows], {k: v.detach().numpy() for k, v in network.state_dict().items()})
            np.testing.assert_allclose(handle['X'][rows], independent, rtol=3e-5, atol=2e-5)
        report['contexts'][context] = {'guide_group_isolation': True, 'heldout_cells': 2800,
            'all_gene_and_cell_mse_recomputed': True, 'full_raw_cells_scanned': 18400,
            'independent_forward_rows_per_model': 32, 'normalized_mean_baseline_genes_checked_per_model': 32,
            'raw_baseline_genes_checked_per_model': 18233, 'zero_input_output_overlap': True,
            'validation_only_selection': True, 'original_cell_ids_and_counts_preserved': True}
        del raw, source
    write_json(out / 'validation.json', report)
    print('Neural model, metrics, source alignment and reconstruction audit passed.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT)
    args = parser.parse_args()
    torch.set_num_threads(8)
    with threadpool_limits(limits=8):
        audit(args.output)
