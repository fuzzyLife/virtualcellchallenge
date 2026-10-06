"""Independently audit saved 300-gene models, splits, metrics and prediction files."""
import argparse
import json
from pathlib import Path
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.linalg import solve
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits
from build_submission import ROOT, read_spec, sha256, write_json
from train_gene300_models import OUT, load_model, predict


def strings(node):
    """Decode modern AnnData string/categorical encodings without loading X."""
    if isinstance(node, h5py.Dataset):
        return node.asstr()[:]
    encoding = node.attrs.get('encoding-type')
    if encoding == 'nullable-string-array':
        if node['mask'][:].any():
            raise ValueError('Missing cell/gene labels')
        return strings(node['values'])
    if encoding == 'categorical':
        codes = node['codes'][:]
        if (codes < 0).any():
            raise ValueError('Missing categorical labels')
        return strings(node['categories'])[codes]
    raise ValueError(f'Unsupported label encoding {encoding}')


def independent_features(inputs, family):
    # Independently written feature construction, without importing the training helper.
    counts = np.asarray(inputs, dtype=np.float64)
    pieces = []
    if family in ('rank', 'hybrid'):
        ranks = np.vstack([rankdata(row, method='average') for row in counts])
        ranks -= 150.5
        norms = np.sqrt((ranks ** 2).sum(axis=1))
        pieces.append(np.divide(ranks, norms[:, None], out=np.zeros_like(ranks), where=norms[:, None] > 0))
    if family in ('magnitude', 'hybrid'):
        pieces.append(np.log1p(counts))
    pieces.append(np.column_stack((np.log1p(counts.sum(axis=1)), (counts > 0).sum(axis=1) / 300)))
    return np.concatenate(pieces, axis=1)


def audit(out):
    manifest, genes, targets = read_spec(ROOT)
    run = json.loads((out / 'run_manifest.json').read_text(encoding='utf-8'))
    assert set(run['contexts']) == {'A', 'B', 'C'}
    for name, digest in run['source_sha256'].items():
        assert sha256(ROOT / name) == digest, name
    assert sha256(ROOT / 'train_gene300_models.py') == run['code_sha256']
    for name, expected in run['outputs'].items():
        assert (out / name).stat().st_size == expected['bytes']
        assert sha256(out / name) == expected['sha256'], name
    splits = pd.read_csv(out / 'cell_splits.csv')
    metrics = pd.read_csv(out / 'per_gene_test_metrics.csv')
    cells = pd.read_csv(out / 'per_cell_test_metrics.csv')
    summary = pd.read_csv(out / 'test_summary.csv').set_index('context')
    candidates = pd.read_csv(out / 'candidate_scores.csv')
    assert len(splits) == 55200 and len(metrics) == 54699 and len(cells) == 8400
    positions = {g: i for i, g in enumerate(genes)}
    ipos = [positions[g] for g in targets]
    opos = [i for i, g in enumerate(genes) if g not in set(targets)]
    report = {'passed': True, 'contexts': {}, 'all_output_hashes_match': True,
              'all_source_hashes_match': True, 'coefficient_check': '24 output genes per model via independent SPD solve',
              'full_prediction_check': 'all values scanned; 24 cell rows per context replayed for all output genes'}
    for context in 'ABC':
        print(f'{context}: validating split isolation, coefficients and every saved expression value', flush=True)
        obj = ad.read_h5ad(ROOT / f'context_{context}.h5ad')
        assert sha256(ROOT / f'context_{context}.h5ad') == run['contexts'][context]['source_sha256']
        raw = obj.X.astype(np.int32).tocsr()
        raw.sum_duplicates(); raw.eliminate_zeros(); raw.sort_indices()
        panel = raw[:, ipos].toarray().astype(np.float64)
        libraries = np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel()
        split = splits[splits.context == context].sort_values('source_row')
        np.testing.assert_array_equal(split.source_row, np.arange(18400))
        np.testing.assert_array_equal(split.cell_id, obj.obs_names)
        np.testing.assert_array_equal(split.ntc_id, obj.obs.ntc_id.astype(str))
        np.testing.assert_array_equal(split.input_300_total, panel.sum(axis=1))
        assert split.groupby('ntc_id')['split'].nunique().max() == 1
        assert split.groupby('split').size().to_dict() == {'train': 12800, 'validation': 2800, 'test': 2800}
        assert split.groupby('split').ntc_id.nunique().to_dict() == {'train': 32, 'validation': 7, 'test': 7}
        train = split['split'].to_numpy() == 'train'
        test = split['split'].to_numpy() == 'test'
        held_rows = np.flatnonzero(test)
        selected = json.loads((out / f'selected_{context}.json').read_text(encoding='utf-8'))
        candidate = candidates[candidates.context == context].sort_values('validation_expressed_logcpm_mse').iloc[0]
        assert selected['family'] == candidate.family and selected['alpha'] == candidate.alpha
        gene_table = metrics[metrics.context == context]
        np.testing.assert_array_equal(gene_table.gene_name, genes[opos])
        chosen_genes = np.linspace(0, len(opos) - 1, 24, dtype=int)
        for prefix, fit_mask, scope in [('evaluation_model', train, 'train_only_evaluation'),
                                       ('model', np.ones(len(panel), dtype=bool), 'all_controls_production')]:
            model = load_model(out / f'{prefix}_{context}.npz')
            assert str(model['context']) == context and str(model['scope']) == scope
            assert str(model['family']) == selected['family'] and float(model['alpha']) == selected['alpha']
            np.testing.assert_array_equal(model['input_genes'], targets)
            np.testing.assert_array_equal(model['output_genes'], genes[opos])
            np.testing.assert_array_equal(model['fit_row_indices'], np.flatnonzero(fit_mask))
            f = independent_features(panel[fit_mask], str(model['family']))
            mean, std = f.mean(axis=0), f.std(axis=0)
            std[std < 1e-12] = 1
            np.testing.assert_allclose(model['feature_mean'], mean, atol=1e-12)
            np.testing.assert_allclose(model['feature_scale'], std, atol=1e-12)
            x = (f - mean) / std
            y = np.log1p(raw[fit_mask][:, np.asarray(opos)[chosen_genes]].toarray() *
                         (1e6 / libraries[fit_mask])[:, None]).astype(np.float32).astype(np.float64)
            intercept = y.mean(axis=0)
            gram = x.T @ x / len(x) + float(model['alpha']) * np.eye(x.shape[1])
            coefficient = solve(gram, x.T @ (y - intercept) / len(x), assume_a='pos')
            np.testing.assert_allclose(model['intercept'][chosen_genes], intercept, atol=1e-10)
            np.testing.assert_allclose(model['coefficient'][:, chosen_genes], coefficient, rtol=1e-5, atol=2e-7)
            fitted = x @ model['coefficient'][:, chosen_genes] + intercept
            np.testing.assert_allclose(model['cpm_smearing'][chosen_genes], np.exp(y - fitted).mean(axis=0), rtol=1e-10)
            dfeatures = independent_features(panel[fit_mask], 'magnitude')
            dmean, dscale = dfeatures.mean(axis=0), dfeatures.std(axis=0)
            dscale[dscale < 1e-12] = 1
            np.testing.assert_allclose(model['depth_feature_mean'], dmean)
            np.testing.assert_allclose(model['depth_feature_scale'], dscale)
            dx = (dfeatures - dmean) / dscale
            ratio = np.log1p(libraries[fit_mask]) - np.log1p(panel[fit_mask].sum(axis=1))
            dc = solve(dx.T @ dx / len(dx) + float(model['depth_alpha']) * np.eye(dx.shape[1]),
                       dx.T @ (ratio - ratio.mean()) / len(dx), assume_a='pos')
            np.testing.assert_allclose(model['depth_coefficient'], dc, rtol=1e-7, atol=1e-9)
            np.testing.assert_allclose(model['depth_intercept'][0], ratio.mean(), atol=1e-12)
            residual = ratio - dx @ model['depth_coefficient'] - ratio.mean()
            np.testing.assert_allclose(model['depth_smearing'], np.exp(residual).mean(), rtol=1e-10)
        evaluation = load_model(out / f'evaluation_model_{context}.npz')
        production = load_model(out / f'model_{context}.npz')
        with h5py.File(out / f'heldout_predictions_{context}.h5ad', 'r') as handle:
            assert handle['X'].shape == (2800, 18233)
            np.testing.assert_array_equal(strings(handle['var/gene_name']), genes[opos])
            np.testing.assert_array_equal(strings(handle['obs/_index']), obj.obs_names[test])
            np.testing.assert_array_equal(handle['obsm/measured_300_raw_counts'][:], panel[test])
            saved_depth = handle['obs/predicted_library_counts'][:]
            cell_table = cells[cells.context == context]
            np.testing.assert_array_equal(cell_table.cell_id, obj.obs_names[test])
            np.testing.assert_allclose(cell_table.predicted_library_counts, saved_depth)
            logmse, rawmse, baselogmse, baserawmse, sizesemse, expressed = [], [], [], [], [], []
            for start in range(0, len(opos), 256):
                stop = min(start + 256, len(opos))
                p = handle['X'][:, start:stop].astype(np.float64)
                pl = handle['layers/predicted_log1p_CPM'][:, start:stop].astype(np.float64)
                assert np.isfinite(p).all() and (p >= 0).all() and np.isfinite(pl).all()
                truth = raw[test][:, opos[start:stop]].toarray().astype(np.float64)
                yl = np.log1p(truth * (1e6 / libraries[test])[:, None]).astype(np.float32).astype(np.float64)
                trainraw = raw[train][:, opos[start:stop]].toarray().astype(np.float64)
                traincpm = trainraw * (1e6 / libraries[train])[:, None]
                trainlogmean = np.log1p(traincpm).astype(np.float32).mean(axis=0, dtype=np.float64)
                expressed.extend(traincpm.mean(axis=0) >= 5)
                baselogmse.extend(((trainlogmean - yl) ** 2).mean(axis=0))
                baserawmse.extend(((trainraw.mean(axis=0) - truth) ** 2).mean(axis=0))
                sizebase = panel[test].sum(axis=1)[:, None] * (trainraw.sum(axis=0) / panel[train].sum())
                sizesemse.extend(((sizebase - truth) ** 2).mean(axis=0))
                rawmse.extend(((p - truth) ** 2).mean(axis=0))
                logmse.extend(((pl - yl) ** 2).mean(axis=0))
            np.testing.assert_allclose(gene_table.test_raw_mse, rawmse, rtol=1e-9, atol=1e-10)
            np.testing.assert_allclose(gene_table.test_logcpm_mse, logmse, rtol=1e-9, atol=1e-10)
            np.testing.assert_allclose(gene_table.test_mean_baseline_logcpm_mse, baselogmse, rtol=1e-9, atol=1e-10)
            np.testing.assert_allclose(gene_table.test_mean_baseline_raw_mse, baserawmse, rtol=1e-8, atol=1e-9)
            np.testing.assert_allclose(gene_table.test_depth_baseline_raw_mse, sizesemse, rtol=1e-9, atol=1e-10)
            np.testing.assert_array_equal(gene_table.expressed_train_5cpm, expressed)
            expressed = np.asarray(expressed)
            expected_skill = 1 - np.sum(np.asarray(logmse)[expressed]) / np.sum(np.asarray(baselogmse)[expressed])
            np.testing.assert_allclose(summary.loc[context, 'test_expressed_logcpm_skill'], expected_skill, atol=1e-12)
            raw_skill = 1 - np.sum(rawmse) / np.sum(sizesemse)
            np.testing.assert_allclose(summary.loc[context, 'test_raw_skill_vs_depth_baseline'], raw_skill, atol=1e-12)
            # Replay every held-out cell, all genes, without reading unobserved values as inputs.
            for start in range(0, len(held_rows), 128):
                stop = min(start + 128, len(held_rows))
                p, pl, d = predict(evaluation, panel[held_rows[start:stop]])
                np.testing.assert_allclose(handle['X'][start:stop], p, rtol=2e-6, atol=1e-6)
                np.testing.assert_allclose(handle['layers/predicted_log1p_CPM'][start:stop], pl, rtol=2e-6, atol=1e-6)
                np.testing.assert_allclose(saved_depth[start:stop], d, rtol=1e-10)
        with h5py.File(out / f'reconstructed_{context}.h5ad', 'r') as handle:
            assert handle['X'].shape == (18400, 18233)
            np.testing.assert_array_equal(strings(handle['var/gene_name']), genes[opos])
            np.testing.assert_array_equal(strings(handle['obs/_index']), obj.obs_names)
            np.testing.assert_array_equal(handle['obsm/measured_300_raw_counts'][:], panel)
            assert set(strings(handle['obs/prediction_scope/categories'])) == {'in_sample_all_controls_refit'}
            for start in range(0, 18400, 256):
                values = handle['X'][start:start + 256]
                assert np.isfinite(values).all() and (values >= 0).all()
                np.testing.assert_allclose(values.sum(axis=1, dtype=np.float64) + panel[start:start + 256].sum(axis=1),
                    handle['obs/predicted_library_counts'][start:start + 256], rtol=1e-6)
            rows = np.linspace(0, 18399, 24, dtype=int)
            p, _, d = predict(production, panel[rows])
            np.testing.assert_allclose(handle['X'][rows], p, rtol=2e-6, atol=1e-6)
        report['contexts'][context] = {'split_isolated_by_guide': True, 'train_cells': 12800,
            'validation_cells': 2800, 'test_cells': 2800, 'full_cells': 18400,
            'input_values_checked': int(panel.size), 'output_values_scanned': int(21200 * 18233),
            'heldout_cells_all_genes_replayed': 2800, 'full_fit_cells_all_genes_replayed': 24,
            'gene_coefficients_independent_solve': 24 * 2, 'depth_coefficients_independent_solve': True,
            'all_heldout_gene_raw_and_logcpm_mse_recomputed': True}
        del obj, raw
    write_json(out / 'validation.json', report)
    print('All saved models, split isolation, held-out MSE values and reconstruction scans passed.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT)
    args = parser.parse_args()
    with threadpool_limits(limits=8):
        audit(args.output)
