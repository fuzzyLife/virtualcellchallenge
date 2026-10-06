"""Separate control-cell expression models: 300 inputs -> 18,233 other genes.
No full-library totals or unobserved genes are used as inference inputs.
"""
from __future__ import annotations
import argparse
import gc
from importlib.metadata import version
from pathlib import Path
import time
import anndata as ad
import numpy as np
import pandas as pd
from scipy.linalg import eigh
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits
from build_submission import ROOT, load_control, read_spec, sha256, write_json

OUT = ROOT / 'analysis/gene300_models'
ALPHAS = (.01, .1, 1., 10.)
FAMILIES = ('rank', 'magnitude', 'hybrid')
SEED = 20261006
BLOCK = 256


def log(message):
    print(time.strftime('%H:%M:%S') + ' ' + message, flush=True)


def split_guides(guides, seed=SEED):
    names = np.unique(np.asarray(guides).astype(str))
    if len(names) != 46:
        raise ValueError('Expected 46 control guides')
    names = np.random.default_rng(seed).permutation(names)
    mapping = {g: role for role, group in zip(('train', 'validation', 'test'),
               (names[:32], names[32:39], names[39:])) for g in group}
    return np.array([mapping[str(g)] for g in guides])


def features(counts, family):
    """Zeros participate in tied ranks; no pseudocount is used for ranking."""
    x = np.asarray(counts, dtype=np.float64)
    if x.ndim != 2 or x.shape[1] != 300 or not np.isfinite(x).all() or (x < 0).any():
        raise ValueError('Expected finite nonnegative counts with exactly 300 columns')
    summaries = np.column_stack((np.log1p(x.sum(axis=1)), (x > 0).mean(axis=1)))
    pieces = []
    if family in ('rank', 'hybrid'):
        ranks = rankdata(x, axis=1, method='average')
        ranks -= ranks.mean(axis=1, keepdims=True)
        norm = np.linalg.norm(ranks, axis=1, keepdims=True)
        pieces.append(np.divide(ranks, norm, out=np.zeros_like(ranks), where=norm > 0))
    if family in ('magnitude', 'hybrid'):
        pieces.append(np.log1p(x))
    if family not in FAMILIES:
        raise ValueError(f'Unknown feature family {family}')
    return np.column_stack((*pieces, summaries))


def standardize(f, mean=None, scale=None):
    if mean is None:
        mean, scale = f.mean(axis=0), f.std(axis=0)
        scale[scale < 1e-12] = 1.
    return (f - mean) / scale, mean, scale


def ridge_basis(x):
    values, vectors = eigh(x.T @ x / len(x))
    return np.maximum(values, 0), vectors, x @ vectors


def ridge_fit(x, y, alpha):
    mean = np.asarray(y, dtype=np.float64).mean(axis=0)
    values, vectors, z = ridge_basis(x)
    cross = z.T @ (y - mean) / len(x)
    return vectors @ (cross / (values[:, None] + alpha)), mean


def validation_sse(b, cross, gram, centered_sse):
    return np.maximum(0., centered_sse - 2 * np.sum(b * cross, axis=0)
                      + np.sum(b * (gram @ b), axis=0))


def fit_candidates(inputs, labels, train, validation, expressed):
    """Evaluate exact ridge validation loss using sufficient statistics."""
    ymean = labels[train].mean(axis=0, dtype=np.float64)
    records = []
    for family in FAMILIES:
        log(f'  selecting {family} features')
        f = features(inputs, family)
        xt, mean, scale = standardize(f[train])
        xv, _, _ = standardize(f[validation], mean, scale)
        values, vectors, zt = ridge_basis(xt)
        zv = xv @ vectors
        gram = zv.T @ zv
        sums = np.zeros((len(ALPHAS), labels.shape[1]))
        baseline = np.empty(labels.shape[1])
        for start in range(0, labels.shape[1], BLOCK):
            stop = min(start + BLOCK, labels.shape[1])
            yt = labels[train, start:stop].astype(np.float64) - ymean[start:stop]
            yv = labels[validation, start:stop].astype(np.float64) - ymean[start:stop]
            cross = zt.T @ yt / len(xt)
            valcross = zv.T @ yv
            base = np.sum(yv * yv, axis=0)
            baseline[start:stop] = base
            for a, alpha in enumerate(ALPHAS):
                b = cross / (values[:, None] + alpha)
                sums[a, start:stop] = validation_sse(b, valcross, gram, base)
        for alpha, sse in zip(ALPHAS, sums):
            records.append({'family': family, 'alpha': alpha, 'features': xt.shape[1],
                'validation_expressed_logcpm_mse': float(sse[expressed].mean() / len(xv)),
                'validation_all_logcpm_mse': float(sse.mean() / len(xv)),
                'validation_expressed_skill': float(1 - sse[expressed].sum() / baseline[expressed].sum())})
    return min(records, key=lambda r: r['validation_expressed_logcpm_mse']), records


def choose_depth(inputs, libraries, train, validation):
    f = features(inputs, 'magnitude')
    x, mean, scale = standardize(f[train])
    xv, _, _ = standardize(f[validation], mean, scale)
    y = np.log1p(libraries) - np.log1p(inputs.sum(axis=1))
    records = []
    for alpha in ALPHAS:
        coef, intercept = ridge_fit(x, y[train, None], alpha)
        mse = float(np.mean((xv @ coef[:, 0] + intercept[0] - y[validation]) ** 2))
        records.append({'alpha': alpha, 'validation_log_depth_ratio_mse': mse})
    return min(records, key=lambda r: r['validation_log_depth_ratio_mse'])['alpha'], records


def fit_model(inputs, labels, libraries, indices, selected, input_genes, output_genes,
              context, scope, depth_alpha):
    f = features(inputs[indices], selected['family'])
    x, mean, scale = standardize(f)
    values, vectors, z = ridge_basis(x)
    ymean = labels[indices].mean(axis=0, dtype=np.float64)
    coef = np.empty((x.shape[1], labels.shape[1]), dtype=np.float32)
    smear, fallback = np.empty(labels.shape[1]), np.empty(labels.shape[1])
    for start in range(0, labels.shape[1], BLOCK):
        stop = min(start + BLOCK, labels.shape[1])
        y = labels[indices, start:stop].astype(np.float64)
        cross = z.T @ (y - ymean[start:stop]) / len(x)
        b = vectors @ (cross / (values[:, None] + selected['alpha']))
        coef[:, start:stop] = b
        fitted = x @ coef[:, start:stop] + ymean[start:stop]
        smear[start:stop] = np.mean(np.exp(y - fitted), axis=0)
        fallback[start:stop] = np.mean(np.expm1(y), axis=0)
    dx, dmean, dscale = standardize(features(inputs[indices], 'magnitude'))
    ratio = np.log1p(libraries[indices]) - np.log1p(inputs[indices].sum(axis=1))
    dc, intercept = ridge_fit(dx, ratio[:, None], depth_alpha)
    resid = ratio - (dx @ dc).ravel() - intercept[0]
    return {'input_genes': np.asarray(input_genes).astype(str), 'output_genes': np.asarray(output_genes).astype(str),
        'context': np.array(context), 'scope': np.array(scope), 'family': np.array(selected['family']),
        'alpha': np.array(selected['alpha']), 'feature_mean': mean, 'feature_scale': scale,
        'coefficient': coef, 'intercept': ymean, 'cpm_smearing': smear, 'fallback_cpm': fallback,
        'depth_alpha': np.array(depth_alpha), 'depth_feature_mean': dmean, 'depth_feature_scale': dscale,
        'depth_coefficient': dc.ravel(), 'depth_intercept': intercept,
        'depth_smearing': np.array(np.mean(np.exp(resid))),
        'depth_ratio_bounds': np.array([ratio.min() - 1, ratio.max() + 1]),
        'fit_row_indices': np.flatnonzero(indices), 'schema': np.array(1)}


def load_model(path):
    with np.load(path, allow_pickle=False) as archive:
        model = {k: archive[k] for k in archive.files}
    if int(model['schema']) != 1:
        raise ValueError('Unsupported model schema')
    return model


def predict(model, inputs):
    """Uses only 300 measured genes, including for inference of total counts."""
    x, _, _ = standardize(features(inputs, str(model['family'])), model['feature_mean'], model['feature_scale'])
    latent = np.clip(x @ model['coefficient'] + model['intercept'], 0, np.log1p(1e6))
    cpm = np.clip(np.exp(latent) * model['cpm_smearing'] - 1., 0., 1e6)
    dx, _, _ = standardize(features(inputs, 'magnitude'), model['depth_feature_mean'], model['depth_feature_scale'])
    ratio = np.clip(dx @ model['depth_coefficient'] + model['depth_intercept'][0], *model['depth_ratio_bounds'])
    measured = inputs.sum(axis=1, dtype=np.float64)
    depth = np.maximum(measured, np.exp(np.log1p(measured) + ratio) * model['depth_smearing'] - 1.)
    total = cpm.sum(axis=1)
    empty = total <= 0
    if empty.any():
        cpm[empty] = model['fallback_cpm']
        total = cpm.sum(axis=1)
    proportions = np.divide(cpm, total[:, None], out=np.zeros_like(cpm), where=total[:, None] > 0)
    counts = proportions * (depth - measured)[:, None]
    return counts.astype(np.float32), latent.astype(np.float32), depth


def correlation(a, b, axis=0):
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    a, b = a - a.mean(axis=axis, keepdims=True), b - b.mean(axis=axis, keepdims=True)
    norm = np.sqrt(np.sum(a * a, axis=axis) * np.sum(b * b, axis=axis))
    return np.divide(np.sum(a * b, axis=axis), norm, out=np.full(norm.shape, np.nan), where=norm > 0)


def spearman(a, b, axis=0):
    # Measured zeros remain valid; ties use average ranks; constant vectors give NaN.
    return correlation(rankdata(a, method='average', axis=axis), rankdata(b, method='average', axis=axis), axis)


def evaluate(model, inputs, raw_outputs, labels, libraries, roles, obs, expressed):
    train, test = roles == 'train', roles == 'test'
    indices = np.flatnonzero(test)
    predicted = np.empty((len(indices), labels.shape[1]), dtype=np.float32)
    latent, depth = np.empty_like(predicted), np.empty(len(indices))
    for start in range(0, len(indices), 128):
        stop = min(start + 128, len(indices))
        predicted[start:stop], latent[start:stop], depth[start:stop] = predict(model, inputs[indices[start:stop]])
    logmean = labels[train].mean(axis=0, dtype=np.float64)
    rawmean = np.asarray(raw_outputs[train].mean(axis=0)).ravel()
    size_coef = np.asarray(raw_outputs[train].sum(axis=0)).ravel() / inputs[train].sum()
    records = []
    raw_cell, size_cell, log_cell, logbase_cell = [np.zeros(len(indices)) for _ in range(4)]
    for start in range(0, labels.shape[1], BLOCK):
        stop = min(start + BLOCK, labels.shape[1])
        truth = raw_outputs[test, start:stop].toarray().astype(np.float64)
        yl = labels[test, start:stop].astype(np.float64)
        p, pl = predicted[:, start:stop].astype(np.float64), latent[:, start:stop].astype(np.float64)
        base_size = inputs[test].sum(axis=1)[:, None] * size_coef[None, start:stop]
        se, se_size = (p - truth) ** 2, (base_size - truth) ** 2
        se_log, se_logbase = (pl - yl) ** 2, (logmean[start:stop] - yl) ** 2
        raw_cell += se.sum(axis=1); size_cell += se_size.sum(axis=1)
        log_cell += se_log.sum(axis=1); logbase_cell += se_logbase.sum(axis=1)
        lm, lb = se_log.mean(axis=0), se_logbase.mean(axis=0)
        records.append(pd.DataFrame({'context': str(model['context']), 'gene_name': model['output_genes'][start:stop],
            'expressed_train_5cpm': expressed[start:stop], 'train_mean_raw': rawmean[start:stop],
            'train_mean_cpm': model['fallback_cpm'][start:stop], 'test_raw_mse': se.mean(axis=0),
            'test_mean_baseline_raw_mse': ((truth - rawmean[start:stop]) ** 2).mean(axis=0),
            'test_depth_baseline_raw_mse': se_size.mean(axis=0), 'test_logcpm_mse': lm,
            'test_mean_baseline_logcpm_mse': lb,
            'test_logcpm_skill': np.divide(lb - lm, lb, out=np.full(len(lm), np.nan), where=lb > 0),
            'test_raw_pearson': correlation(p, truth), 'test_normalized_spearman': spearman(pl, yl)}))
    genes = pd.concat(records, ignore_index=True)
    cell_corr, base_corr, size_corr = [], [], []
    for start in range(0, len(indices), 128):
        stop = min(start + 128, len(indices))
        truth = raw_outputs[indices[start:stop]].toarray()
        cell_corr.extend(spearman(predicted[start:stop], truth, axis=1))
        base_corr.extend(spearman(np.broadcast_to(rawmean, truth.shape), truth, axis=1))
        size_corr.extend(spearman(inputs[indices[start:stop]].sum(axis=1)[:, None] * size_coef, truth, axis=1))
    cells = obs.iloc[indices].copy()
    cells['source_row'] = indices
    cells['predicted_library_counts'] = depth
    cells['actual_library_counts_evaluation_only'] = libraries[test]
    for name, values in [('raw_mse', raw_cell), ('depth_baseline_raw_mse', size_cell),
                         ('logcpm_mse', log_cell), ('mean_baseline_logcpm_mse', logbase_cell)]:
        cells[name] = values / labels.shape[1]
    cells['raw_profile_spearman'] = cell_corr
    cells['mean_baseline_profile_spearman'] = base_corr
    cells['depth_baseline_profile_spearman'] = size_corr
    relevant = genes[genes.expressed_train_5cpm]
    summary = {'context': str(model['context']), 'selected_family': str(model['family']), 'alpha': float(model['alpha']),
        'test_cells': len(indices), 'output_genes': len(genes), 'expressed_genes': int(expressed.sum()),
        'test_expressed_logcpm_mse': float(relevant.test_logcpm_mse.mean()),
        'test_expressed_logcpm_skill': float(1 - relevant.test_logcpm_mse.sum() / relevant.test_mean_baseline_logcpm_mse.sum()),
        'test_all_logcpm_skill': float(1 - genes.test_logcpm_mse.sum() / genes.test_mean_baseline_logcpm_mse.sum()),
        'test_raw_mse': float(genes.test_raw_mse.mean()),
        'test_raw_skill_vs_mean': float(1 - genes.test_raw_mse.sum() / genes.test_mean_baseline_raw_mse.sum()),
        'test_raw_skill_vs_depth_baseline': float(1 - genes.test_raw_mse.sum() / genes.test_depth_baseline_raw_mse.sum()),
        'median_expressed_gene_normalized_spearman': float(relevant.test_normalized_spearman.median()),
        'median_cell_profile_spearman': float(np.nanmedian(cell_corr)),
        'median_mean_baseline_cell_profile_spearman': float(np.nanmedian(base_corr)),
        'median_library_absolute_relative_error': float(np.median(np.abs(depth - libraries[test]) / libraries[test]))}
    return summary, genes, cells, predicted, latent, depth


def write_predictions(path, model, inputs, obs, scope, predictions=None, latent=None, depth=None):
    if predictions is None:
        predictions = np.empty((len(inputs), len(model['output_genes'])), dtype=np.float32)
        depth = np.empty(len(inputs))
        for start in range(0, len(inputs), 128):
            stop = min(start + 128, len(inputs))
            predictions[start:stop], _, depth[start:stop] = predict(model, inputs[start:stop])
    metadata = obs.copy()
    metadata['prediction_scope'], metadata['predicted_library_counts'] = scope, depth
    obj = ad.AnnData(predictions, obs=metadata, var=pd.DataFrame(index=pd.Index(model['output_genes'], name='gene_name')))
    obj.obsm['measured_300_raw_counts'] = np.asarray(inputs, dtype=np.int32)
    obj.uns['input_genes'] = model['input_genes']
    obj.uns['units'] = 'nonnegative floating expected raw counts'
    obj.uns['model_family'], obj.uns['prediction_scope'] = str(model['family']), scope
    if latent is not None:
        obj.layers['predicted_log1p_CPM'] = latent
    obj.write_h5ad(path, compression='lzf')


def run(args):
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    manifest, genes, targets = read_spec(ROOT)
    from analyze_target_archive import read_archive
    aframe, ainputs, aaudit = read_archive(ROOT / 'vcc_2026_300_target_matrix.zip', targets)
    positions = {g: i for i, g in enumerate(genes)}
    ipos = np.array([positions[g] for g in targets])
    opos = np.array([i for i, g in enumerate(genes) if g not in set(targets)])
    output_genes = genes[opos]
    sources = ['manifest.json', 'gene_names.csv', 'pert_counts.csv', 'vcc_2026_300_target_matrix.zip']
    provenance = {'seed': SEED, 'alphas': list(ALPHAS), 'families': list(FAMILIES),
        'primary_selection': 'validation log1p CPM MSE, training mean CPM >=5',
        'split_unit': 'NTC guide ID; 32/7/7 train/validation/test',
        'source_sha256': {p: sha256(ROOT / p) for p in sources}, 'archive': aaudit, 'contexts': {},
        'versions': {p: version(p) for p in ('numpy', 'scipy', 'pandas', 'anndata', 'h5py', 'threadpoolctl')},
        'code_sha256': sha256(Path(__file__)), 'outputs': {}}
    summaries, gtab, ctab, splits, candidates, dcandidates = [], [], [], [], [], []
    for context in args.contexts:
        log(f'{context}: loading and validating every control count')
        raw, guides, libraries, _, control_summary, _ = load_control(ROOT, context, genes, manifest)
        source = ad.read_h5ad(ROOT / f'context_{context}.h5ad', backed='r')
        obs = source.obs.copy(); source.file.close()
        obs['cell_id'], obs['ntc_id'], obs['context'] = obs.index.astype(str), np.asarray(guides).astype(str), context
        if not obs.index.is_unique:
            raise ValueError('Duplicate source cell IDs')
        inputs = raw[:, ipos].toarray().astype(np.float64)
        if not np.array_equal(inputs, ainputs[aframe.context.to_numpy() == context]):
            raise ValueError(f'{context}: ZIP inputs disagree with original H5AD')
        roles = split_guides(guides)
        train, val, test = roles == 'train', roles == 'validation', roles == 'test'
        split = obs[['cell_id', 'context', 'ntc_id']].copy()
        split['source_row'], split['split'], split['input_300_total'] = np.arange(raw.shape[0]), roles, inputs.sum(axis=1)
        splits.append(split)
        other = raw[:, opos].tocsr()
        labels, mean_cpm = np.empty(other.shape, dtype=np.float32), np.empty(other.shape[1])
        for start in range(0, other.shape[1], BLOCK):
            stop = min(start + BLOCK, other.shape[1])
            cpm = other[:, start:stop].toarray().astype(np.float64) * (1e6 / libraries)[:, None]
            labels[:, start:stop] = np.log1p(cpm)
            mean_cpm[start:stop] = cpm[train].mean(axis=0)
        expressed = mean_cpm >= 5
        selected, scores = fit_candidates(inputs, labels, train, val, expressed)
        da, dscores = choose_depth(inputs, libraries, train, val)
        candidates.extend(dict(context=context, **r) for r in scores)
        dcandidates.extend(dict(context=context, **r) for r in dscores)
        # All settings are locked before test labels are evaluated.
        write_json(out / f'selected_{context}.json', dict(context=context, **selected, depth_alpha=da))
        log(f'{context}: selected {selected["family"]}, alpha={selected["alpha"]}; fitting evaluation model')
        model = fit_model(inputs, labels, libraries, train, selected, targets, output_genes, context, 'train_only_evaluation', da)
        np.savez_compressed(out / f'evaluation_model_{context}.npz', **model)
        summary, gs, cs, p, lp, depth = evaluate(model, inputs, other, labels, libraries, roles, obs, expressed)
        summaries.append(summary); gtab.append(gs); ctab.append(cs)
        log(f'{context}: test normalized skill={summary["test_expressed_logcpm_skill"]:.4f}, '
            f'raw skill vs input-depth baseline={summary["test_raw_skill_vs_depth_baseline"]:.4f}')
        write_predictions(out / f'heldout_predictions_{context}.h5ad', model, inputs[test], obs.iloc[np.flatnonzero(test)],
                          'unseen_test_guides_train_only_model', p, lp, depth)
        del p, lp, model
        log(f'{context}: refitting on all 18,400 controls')
        production = fit_model(inputs, labels, libraries, np.ones(len(inputs), dtype=bool), selected, targets,
                               output_genes, context, 'all_controls_production', da)
        np.savez_compressed(out / f'model_{context}.npz', **production)
        log(f'{context}: saving fitted reconstructions for every cell')
        write_predictions(out / f'reconstructed_{context}.h5ad', production, inputs, obs, 'in_sample_all_controls_refit')
        provenance['contexts'][context] = {'source_sha256': sha256(ROOT / f'context_{context}.h5ad'),
            'control_summary': control_summary, 'archive_values_checked': int(inputs.size), 'archive_mismatches': 0,
            'split_cells': {r: int((roles == r).sum()) for r in ('train', 'validation', 'test')},
            'selected': selected, 'depth_alpha': da, 'summary': summary}
        for prefix, suffix in (('evaluation_model', '.npz'), ('model', '.npz'),
                               ('heldout_predictions', '.h5ad'), ('reconstructed', '.h5ad')):
            name = f'{prefix}_{context}{suffix}'
            provenance['outputs'][name] = {'bytes': (out / name).stat().st_size, 'sha256': sha256(out / name)}
        del raw, other, labels, inputs, production
        gc.collect()
        for name, frames in [('candidate_scores.csv', candidates), ('depth_candidates.csv', dcandidates), ('test_summary.csv', summaries)]:
            pd.DataFrame(frames).to_csv(out / name, index=False)
        pd.concat(gtab, ignore_index=True).to_csv(out / 'per_gene_test_metrics.csv', index=False)
        pd.concat(ctab, ignore_index=True).to_csv(out / 'per_cell_test_metrics.csv', index=False)
        pd.concat(splits, ignore_index=True).to_csv(out / 'cell_splits.csv', index=False)
        write_json(out / 'run_manifest.json', provenance)
    log('Finished separate context models and reconstructions')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT)
    parser.add_argument('--contexts', nargs='+', choices=['A', 'B', 'C'], default=['A', 'B', 'C'])
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()
    with threadpool_limits(limits=args.threads):
        run(args)
