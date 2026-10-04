"""Raw-count, pairwise-complete Spearman correlations on the 300 VCC targets.

No pseudocounts, log transforms, normalization, zero imputation or positive-only
filtering. Measured zero is valid. A missing source symbol is NaN. Each cell is
compared with each of the five previously prepared non-targeting profiles.
"""
from pathlib import Path
from datetime import datetime, timezone
from importlib.metadata import version
import argparse
import hashlib
import json
import platform

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parent
REFERENCES = ['K562', 'HCT116', 'HEK293T', 'H1', 'CD4']
OUT = ROOT / 'analysis/target_rank_correlations'


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b''):
            result.update(block)
    return result.hexdigest()


def pairwise_spearman(left, right, min_pairs=2):
    """Remove nonfinite pairs BEFORE ranking; return rho, n, and status."""
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    if left.shape != right.shape or left.ndim != 1:
        raise ValueError('Expected equally sized one-dimensional vectors')
    valid = np.isfinite(left) & np.isfinite(right)
    count = int(valid.sum())
    if count < min_pairs:
        return np.nan, count, 'insufficient_pairs'
    a = rankdata(left[valid], method='average')
    b = rankdata(right[valid], method='average')
    a -= a.mean()
    b -= b.mean()
    denominator = np.sqrt(np.dot(a, a) * np.dot(b, b))
    if denominator == 0:
        return np.nan, count, 'constant_rank_vector'
    return float(np.clip(np.dot(a, b) / denominator, -1, 1)), count, 'ok'


def pairwise_rank_rows(cells, reference, min_pairs=2):
    """Batch complete rows; re-mask/re-rank individually when a cell has gaps."""
    cells = np.asarray(cells, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    if cells.ndim != 2 or reference.shape != (cells.shape[1],):
        raise ValueError('Cell and reference axes differ')
    mask = np.isfinite(reference)
    subset = cells[:, mask]
    counts = np.isfinite(subset).sum(axis=1)
    rho = np.full(len(cells), np.nan)
    statuses = np.full(len(cells), 'insufficient_pairs', dtype=object)
    if int(mask.sum()) >= min_pairs:
        complete = np.isfinite(subset).all(axis=1)
        rows = np.flatnonzero(complete)
        if len(rows):
            a = rankdata(subset[rows], axis=1, method='average')
            b = rankdata(reference[mask], method='average')
            a -= a.mean(axis=1, keepdims=True)
            b -= b.mean()
            denominator = np.sqrt((a * a).sum(axis=1) * np.dot(b, b))
            valid = denominator > 0
            rho[rows[valid]] = np.clip((a[valid] @ b) / denominator[valid], -1, 1)
            statuses[rows] = np.where(valid, 'ok', 'constant_rank_vector')
        for row in np.flatnonzero(~complete):
            rho[row], counts[row], statuses[row] = pairwise_spearman(cells[row], reference, min_pairs)
    return rho, counts, statuses


def table(frame):
    lines = ['| ' + ' | '.join(frame.columns) + ' |', '| ' + ' | '.join(['---'] * len(frame.columns)) + ' |']
    for row in frame.itertuples(index=False, name=None):
        lines.append('| ' + ' | '.join(f'{x:.4f}' if isinstance(x, (float, np.floating)) else str(x) for x in row) + ' |')
    return '\n'.join(lines)


def write_results(summary, report, pooled):
    medians = summary.pivot(index='context', columns='reference', values='median_spearman').reindex(columns=REFERENCES).reset_index()
    coverage = pd.DataFrame([{'reference': name, 'valid_target_genes': report['reference_valid_target_genes'][name]} for name in REFERENCES])
    text = '# Target-gene rank correlation results\n\n'
    text += f"Generated {report['created_utc']}. All {report['competition_cells']:,} competition controls are included. Gene selection is exactly the 300 entries in `pert_counts.csv`, retaining its order.\n\n"
    text += '## Median per-cell Spearman correlations\n\n' + table(medians) + '\n\n'
    text += 'Each correlation uses the original cell counts and the reference pooled control count vector. Only gene pairs with finite values in both vectors are ranked. Average ranks handle ties, including zeros. There is no added constant, log transform, normalization or missing-value imputation.\n\n'
    text += '## Valid target-gene coverage\n\n' + table(coverage) + '\n\n'
    text += 'K562 is measured on a restricted source gene axis. Its absent targets stay missing and are excluded pairwise. Different reference pairs therefore use different gene subsets. Coverage must be considered when comparing coefficients.\n\n'
    text += '## Pooled A/B/C and reference profiles\n\n' + table(pooled.reset_index()) + '\n\n'
    text += 'Pooled correlations rank summed target-gene counts across controls. They are different statistics from medians of single-cell correlations. The corresponding valid-pair counts are in `pooled_profile_valid_gene_counts.csv`. Individual A/B/C cells are not paired with cells in other contexts.\n\n'
    text += '## Files and reproduction\n\nAll scripts are in the root. See [TARGET_RANK_CORRELATIONS.md](TARGET_RANK_CORRELATIONS.md) for exact commands, schemas, source scope and validation. Numeric results and provenance are in `analysis/target_rank_correlations/`. The previous broad-gene Pearson analysis remains available separately.\n\n'
    text += 'These are unperturbed expression associations restricted to the requested target panel. They do not measure the accuracy of CRISPRi response transfer or establish anonymous cell identities.\n'
    (ROOT / 'TARGET_RANK_CORRELATION_RESULTS.md').write_text(text, encoding='utf-8')


def run(batch_size=1024, min_pairs=2):
    OUT.mkdir(parents=True, exist_ok=True)
    spec = json.loads((ROOT / 'manifest.json').read_text())
    genes = pd.read_csv(ROOT / 'gene_names.csv').gene_name.astype(str).to_numpy()
    targets = pd.read_csv(ROOT / 'pert_counts.csv').target_gene.astype(str).to_numpy()
    if len(targets) != 300 or len(set(targets)) != 300:
        raise ValueError('Expected 300 unique target genes')
    gene_index = {gene: i for i, gene in enumerate(genes)}
    positions = np.array([gene_index[gene] for gene in targets])
    profiles = {}
    records = {}
    for name in REFERENCES:
        path = ROOT / f'reference/control_profiles/{name}.csv'
        record = json.loads(path.with_suffix('.json').read_text())
        if digest(path) != record['profile_sha256']:
            raise ValueError(f'{name} profile checksum mismatch')
        source = pd.read_csv(path).set_index('gene_name').control_count_sum
        if not source.index.is_unique:
            raise ValueError('Source symbols must be unique after preparation')
        profiles[name] = source.reindex(targets).to_numpy(dtype=np.float64)
        records[name] = record
    reference_table = pd.DataFrame({'gene_name': targets, **profiles})
    reference_table.to_csv(OUT / 'target_reference_raw_profiles.csv', index=False)
    coverage = pd.DataFrame({'gene_name': targets})
    for name in REFERENCES:
        coverage['valid_' + name] = np.isfinite(profiles[name])
    coverage.to_csv(OUT / 'target_gene_coverage.csv', index=False)
    frames = []
    summaries = []
    guide_summaries = []
    pooled = {}
    inputs = []
    for context in spec['contexts']:
        path = ROOT / f'context_{context}.h5ad'
        obj = ad.read_h5ad(path, backed='r')
        if obj.n_obs != spec['per_context'][context]['control_cells'] or not np.array_equal(obj.var_names.to_numpy(), genes):
            raise ValueError('Competition axes differ')
        if set(obj.obs.target_gene.astype(str)) != {spec['control_label']} or set(obj.obs.context.astype(str)) != {context}:
            raise ValueError('Expected non-targeting controls in the named context')
        cells = obj.n_obs
        scores = np.full((cells, 5), np.nan)
        valid_counts = np.zeros((cells, 5), dtype=np.int64)
        statuses = np.empty((cells, 5), dtype=object)
        libraries = np.zeros(cells)
        target_nonzeros = np.zeros(cells, dtype=np.int64)
        target_finite = np.zeros(cells, dtype=np.int64)
        total_targets = np.zeros(300)
        for start in range(0, cells, batch_size):
            stop = min(start + batch_size, cells)
            x = sparse.csr_matrix(obj.X[start:stop], dtype=np.float64)
            if not np.isfinite(x.data).all() or (x.data < 0).any() or (x.data != np.floor(x.data)).any():
                raise ValueError('Competition raw counts violate the official count specification')
            x.sum_duplicates()
            x.eliminate_zeros()
            libraries[start:stop] = np.asarray(x.sum(axis=1)).ravel()
            raw = x[:, positions].toarray()
            target_nonzeros[start:stop] = np.count_nonzero(raw, axis=1)
            target_finite[start:stop] = np.isfinite(raw).sum(axis=1)
            total_targets += raw.sum(axis=0)
            for j, name in enumerate(REFERENCES):
                scores[start:stop, j], valid_counts[start:stop, j], statuses[start:stop, j] = pairwise_rank_rows(raw, profiles[name], min_pairs)
        frame = pd.DataFrame({'context': context, 'cell_id': obj.obs_names.astype(str), 'ntc_id': obj.obs.ntc_id.astype(str).to_numpy(), 'observed_target_gene': spec['control_label'], 'total_counts': libraries.astype(np.int64), 'finite_target_genes': target_finite, 'nonzero_target_genes': target_nonzeros})
        for j, name in enumerate(REFERENCES):
            frame['spearman_' + name] = scores[:, j]
            frame['valid_gene_pairs_' + name] = valid_counts[:, j]
            frame['status_' + name] = statuses[:, j]
            values = scores[np.isfinite(scores[:, j]), j]
            quantiles = np.quantile(values, [.05, .25, .5, .75, .95]) if len(values) else np.full(5, np.nan)
            summaries.append({'context': context, 'reference': name, 'cells': cells, 'defined_correlations': len(values), 'undefined_correlations': cells-len(values), 'valid_genes_min': int(valid_counts[:, j].min()), 'valid_genes_max': int(valid_counts[:, j].max()), 'mean_spearman': float(values.mean()) if len(values) else np.nan, 'p05_spearman': quantiles[0], 'p25_spearman': quantiles[1], 'median_spearman': quantiles[2], 'p75_spearman': quantiles[3], 'p95_spearman': quantiles[4]})
            for guide, group in frame.groupby('ntc_id', observed=True):
                guide_summaries.append({'context': context, 'ntc_id': guide, 'reference': name, 'cells': len(group), 'median_spearman': group['spearman_' + name].median(), 'defined_correlations': int(group['spearman_' + name].notna().sum())})
        frames.append(frame)
        pooled['context_' + context] = total_targets
        obj.file.close()
        inputs.append({'file': path.name, 'sha256': digest(path), 'cells': cells})
        print(context, 'finished all', cells, 'cells', flush=True)
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(OUT / 'per_cell_target_rank_correlations.csv', index=False)
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT / 'context_reference_summary.csv', index=False)
    pd.DataFrame(guide_summaries).to_csv(OUT / 'guide_reference_summary.csv', index=False)
    pooled.update(profiles)
    names = list(pooled)
    correlation = pd.DataFrame(np.nan, index=names, columns=names)
    counts = pd.DataFrame(0, index=names, columns=names)
    for left in names:
        for right in names:
            correlation.loc[left, right], counts.loc[left, right], _ = pairwise_spearman(pooled[left], pooled[right], min_pairs)
    correlation.to_csv(OUT / 'pooled_profile_rank_correlations.csv', index_label='profile')
    counts.to_csv(OUT / 'pooled_profile_valid_gene_counts.csv', index_label='profile')
    pd.DataFrame({'gene_name': targets, **pooled}).to_csv(OUT / 'pooled_raw_target_profiles.csv', index=False)
    report = {'created_utc': datetime.now(timezone.utc).isoformat(), 'competition_cells': len(combined), 'requested_target_genes': len(targets), 'gene_order': 'pert_counts.csv original order', 'method': 'Pairwise-complete Spearman on original counts, average ranks for ties', 'min_valid_pairs': min_pairs, 'zero_policy': 'Measured zeros retained as valid tied observations', 'missing_policy': 'Unmeasured source symbols and nonfinite pairs omitted before computing ranks; no imputation', 'transforms': [], 'reference_scope': 'Five previously prepared pooled non-targeting profiles; no new reference sampling', 'reference_valid_target_genes': {name: int(np.isfinite(profiles[name]).sum()) for name in REFERENCES}, 'reference_provenance': records, 'competition_inputs': inputs, 'target_list_sha256': digest(ROOT / 'pert_counts.csv'), 'script_sha256': digest(Path(__file__)), 'versions': {name: version(name) for name in ['numpy', 'scipy', 'pandas', 'anndata', 'h5py']}, 'python': platform.python_version(), 'all_expected_rows_present': len(combined) == sum(spec['per_context'][c]['control_cells'] for c in spec['contexts']), 'unique_context_cell_ids': not combined[['context', 'cell_id']].duplicated().any(), 'undefined_per_cell_correlations': {name: int(combined['spearman_' + name].isna().sum()) for name in REFERENCES}, 'outputs': {path.name: {'bytes': path.stat().st_size, 'sha256': digest(path)} for path in sorted(OUT.glob('*.csv'))}}
    (OUT / 'analysis_report.json').write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    write_results(summary, report, correlation.rename_axis('profile'))
    print(summary.pivot(index='context', columns='reference', values='median_spearman').to_string(), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch-size', type=int, default=1024)
    parser.add_argument('--min-valid-pairs', type=int, default=2)
    args = parser.parse_args()
    if args.batch_size < 1 or args.min_valid_pairs < 2:
        parser.error('batch size must be positive and min-valid-pairs at least two')
    run(args.batch_size, args.min_valid_pairs)
