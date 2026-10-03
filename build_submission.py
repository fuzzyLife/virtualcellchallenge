"""Analyse the VCC 2026 bundle and produce a validated upload-ready prediction.vcc.

Run with .competition_env/Scripts/python.exe build_submission.py
The model transfers experimentally measured genome-wide K562 CRISPR responses to each
anonymous context. No challenge perturbation-response labels are available.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import gc
import hashlib
from html.parser import HTMLParser
import importlib.metadata
import json
from pathlib import Path
import re
import subprocess
import sys
import time

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parent
MAX_NNZ = 4_750_000_000
MAX_COUNTS = 1_000_000


def log(message):
    print(time.strftime('%H:%M:%S') + ' ' + message, flush=True)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(8*1024*1024), b''):
            h.update(b)
    return h.hexdigest()


class PageText(HTMLParser):
    def __init__(self):
        super().__init__()
        self.skip = 0
        self.parts = []
    def handle_starttag(self, tag, attrs):
        if tag in ('script', 'style'):
            self.skip += 1
        if tag in ('p', 'div', 'h1', 'h2', 'h3', 'li', 'tr', 'br'):
            self.parts.append('\n')
        if tag == 'a':
            for k, v in attrs:
                if k == 'href' and v and v.startswith('http'):
                    self.parts.append(f' [Link: {v}] ')
    def handle_endtag(self, tag):
        if tag in ('script', 'style'):
            self.skip = max(0, self.skip-1)
    def handle_data(self, data):
        if not self.skip:
            self.parts.append(data)


def inventory(root, analysis):
    """Inspect and hash every original project file; distinguish web assets/data."""
    ignored = {'.competition_env', 'reference', 'analysis', '__pycache__', '.git', '.agents', '.codex'}
    created = {'build_submission.py', 'download_reference.py', 'test_submission.py',
               'requirements-competition.txt', 'COMPETITION_ANALYSIS.md',
               'prediction.h5ad', 'prediction.vcc', 'prediction.h5ad.partial',
               'prediction.vcc.partial', 'RUN_SUBMISSION.ps1'}
    rows = []
    from pypdf import PdfReader
    for path in sorted(root.rglob('*')):
        rel = path.relative_to(root)
        if not path.is_file() or any(x in ignored for x in rel.parts) or path.name in created:
            continue
        if path.suffix == '.h5ad':
            role = 'challenge control counts; fully analysed'
        elif path.name in ('manifest.json', 'gene_names.csv', 'pert_counts.csv'):
            role = 'authoritative panel specification; fully parsed'
        elif path.suffix.lower() == '.pdf':
            reader = PdfReader(path)
            text = '\n\n'.join(p.extract_text() or '' for p in reader.pages)
            (analysis / (path.stem + '.txt')).write_text(text, encoding='utf-8')
            role = 'unrelated bacteriophage genome-generation paper; excluded from modeling'
        elif path.parent == root and path.suffix.lower() == '.html':
            parser = PageText()
            parser.feed(path.read_text(encoding='utf-8', errors='replace'))
            text = re.sub(r'\n\s*\n', '\n', ''.join(parser.parts))
            (analysis / (path.stem + '.txt')).write_text(text, encoding='utf-8')
            role = 'saved competition documentation' if 'Cell.html' not in path.name else '2026 Cell commentary'
        else:
            # Browser JavaScript, styles, SVGs and images are page dependencies,
            # not prediction/training code. Inspect text assets for provenance.
            role = 'saved website asset; not model code or training data'
        row = {'path': rel.as_posix(), 'bytes': path.stat().st_size,
               'sha256': sha256(path), 'role': role}
        if path.suffix.lower() in ('.download', '.css', '.svg', '.html', '.txt') and path.stat().st_size < 10_000_000:
            text = path.read_text(encoding='utf-8', errors='replace')
            row['keyword_hits'] = {x: len(re.findall(x, text, re.I)) for x in
                                   ['perturb', 'submission', 'raw counts', 'context', 'leaderboard']}
        rows.append(row)
    write_json(analysis / 'file_inventory.json', rows)
    return rows


def read_spec(root):
    manifest = json.loads((root / 'manifest.json').read_text())
    genes = pd.read_csv(root / 'gene_names.csv')['gene_name'].astype(str).to_numpy()
    targets = pd.read_csv(root / 'pert_counts.csv')['target_gene'].astype(str).to_numpy()
    assert len(genes) == manifest['n_genes'] and len(set(genes)) == len(genes)
    assert len(targets) == manifest['n_constructs'] and len(set(targets)) == len(targets)
    assert manifest['control_label'] not in targets
    missing = sorted(set(targets)-set(genes))
    if missing:
        raise ValueError(f'Targets missing from gene axis: {missing}')
    return manifest, genes, targets


def quantiles(a):
    return dict(zip(['min', 'p05', 'p25', 'median', 'p75', 'p95', 'max'],
                    [float(x) for x in np.quantile(a, [0, .05, .25, .5, .75, .95, 1])]))


def load_control(root, context, genes, manifest):
    obj = ad.read_h5ad(root / f'context_{context}.h5ad')
    if not np.array_equal(obj.var_names.to_numpy(), genes):
        raise ValueError(f'Gene axis differs in context {context}')
    assert obj.n_obs == manifest['per_context'][context]['control_cells']
    assert set(obj.obs['target_gene'].astype(str)) == {manifest['control_label']}
    assert set(obj.obs['context'].astype(str)) == {context}
    x = obj.X.tocsr()
    for start in range(0, x.nnz, 4_000_000):
        values = x.data[start:start+4_000_000]
        if not (np.isfinite(values).all() and (values >= 0).all() and (values == np.floor(values)).all()):
            raise ValueError(f'{context} is not raw integer counts')
    x = x.astype(np.int32)
    x.sum_duplicates()
    x.eliminate_zeros()
    x.sort_indices()
    libraries = np.asarray(x.sum(axis=1, dtype=np.int64)).ravel()
    if (libraries <= 0).any() or libraries.max() > MAX_COUNTS:
        raise ValueError('Control library size outside supported range')
    guide_names = obj.obs['ntc_id'].astype(str).to_numpy()
    guide_counts = Counter(guide_names)
    assert len(guide_counts) == manifest['per_context'][context]['n_ntc_ids']
    assert set(guide_counts.values()) == {manifest['cells_per_pert']}
    mean = np.asarray(x.mean(axis=0)).ravel()
    sumsquares = np.bincount(x.indices, weights=x.data.astype(np.float64)**2, minlength=len(genes))
    variance = np.maximum(0, sumsquares/x.shape[0]-mean**2)
    detected = np.bincount(x.indices, minlength=len(genes))/x.shape[0]
    # Used only for comparing basal profiles; submission stays in raw counts.
    mean_cpm = np.zeros(len(genes), dtype=np.float64)
    for start in range(0, x.shape[0], 512):
        sub = x[start:start+512].astype(np.float64)
        mean_cpm += np.asarray(sub.multiply((1e6/libraries[start:start+512])[:, None]).sum(axis=0)).ravel()
    mean_cpm /= x.shape[0]
    summary = {'context': context, 'cells': x.shape[0], 'genes': x.shape[1],
               'stored_nonzeros': int(x.nnz), 'sparsity': float(1-x.nnz/np.prod(x.shape)),
               'n_guides': len(guide_counts), 'guide_cell_counts': dict(guide_counts),
               'library_counts': quantiles(libraries), 'genes_detected_per_cell': quantiles(np.diff(x.indptr)),
               'genes_detected_any_cell': int((detected > 0).sum())}
    stats = pd.DataFrame({'context': context, 'gene_name': genes, 'mean_raw_count': mean,
                          'std_raw_count': np.sqrt(variance), 'detection_fraction': detected,
                          'mean_cpm': mean_cpm})
    del obj
    return x, guide_names, libraries, mean_cpm, summary, stats


def guide_stratified_sample(guides, libraries, n, rng):
    """Balance guides and cover each guide's library-size quantiles without replacement."""
    unique = np.unique(guides)
    base, remainder = divmod(n, len(unique))
    extras = set(rng.choice(len(unique), size=remainder, replace=False).tolist())
    selected = []
    for i, guide in enumerate(unique):
        count = base + (i in extras)
        candidates = np.flatnonzero(guides == guide)
        candidates = candidates[np.argsort(libraries[candidates], kind='stable')]
        if count > len(candidates):
            raise ValueError('Too few input cells for stratified sampling')
        if count:
            # Integer edges keep adjacent quantile bins disjoint. Flooring
            # jittered continuous bins can select their shared boundary twice.
            edges = np.floor(np.linspace(0, len(candidates), count+1)).astype(int)
            positions = rng.integers(edges[:-1], edges[1:])
            selected.extend(candidates[positions])
    selected = np.asarray(selected, dtype=np.int64)
    rng.shuffle(selected)
    assert len(selected) == n and len(np.unique(selected)) == n
    return selected


def reference_model(root, genes, targets, analysis):
    """Transfer verified genome-wide K562 pseudobulk responses, with source-construct checks."""
    source = root/'reference'/'K562_gwps_raw_bulk_01.h5ad'
    with h5py.File(source, 'r') as f:
        labels = f['obs/gene_transcript'].asstr()[:]
        cells = f['obs/num_cells_filtered'][:].astype(np.float64)
        fold_expr = f['obs/fold_expr'][:].astype(np.float64)
        core = f['obs/core_control'][:].astype(bool)
        column = f['var/gene_name']
        # The 2022 file uses old HDF5 categorical encoding: codes plus an
        # object reference to the categories dataset, not modern var_names.
        source_symbols = f[column.attrs['categories']].asstr()[:][column[:]]
        matrix = f['X'][:].astype(np.float64)
    if not np.isfinite(matrix).all() or (matrix < 0).any():
        raise ValueError('Public raw pseudobulk reference has invalid means')
    source_targets = np.array([label.split('_', 2)[1] if label.split('_', 1)[0].isdigit()
                               and len(label.split('_', 2)) == 3 else '' for label in labels])
    valid = np.isfinite(cells) & (cells >= 30)
    controls = core & (source_targets == 'non-targeting') & valid
    if not controls.any():
        raise ValueError('No valid non-targeting core controls in the public reference')
    gene_positions = {g:i for i,g in enumerate(genes)}
    shared = np.array(sorted(set(source_symbols)&set(genes)))
    shared_positions = {g:i for i,g in enumerate(shared)}
    dest_positions = np.array([gene_positions[g] for g in shared])
    # Sum duplicate symbols in the source rather than silently selecting a
    # single Ensembl record. The challenge gene axis is unique.
    valid_columns = np.array([g in shared_positions for g in source_symbols])
    mapping = sparse.csr_matrix((np.ones(valid_columns.sum()),
                (np.flatnonzero(valid_columns), [shared_positions[g] for g in source_symbols[valid_columns]])),
                shape=(len(source_symbols), len(shared)))
    matrix = mapping.T.dot(matrix.T).T
    basal_mean = np.average(matrix[controls], axis=0, weights=cells[controls])
    basal = basal_mean*1e6/basal_mean.sum()
    # Only experimentally effective source knockdowns contribute. A 50%
    # source reduction threshold is looser than the challenge's >80% criterion.
    effective = valid & np.isfinite(fold_expr) & (fold_expr <= .5) & (source_targets != 'non-targeting')
    def response(row_indices):
        mean = np.average(matrix[row_indices], axis=0, weights=cells[row_indices])
        cpm = mean*1e6/mean.sum()
        effect = np.log2((cpm+1)/(basal+1))
        n = float(cells[row_indices].sum())
        control_n = float(cells[controls].sum())
        # Suppress small source changes compatible with Poisson count noise.
        # This is a conservative noise approximation, not a DE significance test.
        noise_variance = (1/np.maximum(n*mean,1) +
                          1/np.maximum(control_n*basal_mean,1))/(np.log(2)**2)
        effect = np.sign(effect)*np.maximum(np.abs(effect)-1.96*np.sqrt(noise_variance),0)
        effect = np.clip(effect,-3,3)
        return effect*n/(n+50), n
    # Independent promoter/transcript-targeting constructs for the same source
    # gene provide a replication proxy. They do NOT validate cross-cell-line
    # generalization; the challenge responses remain wholly withheld.
    grouped = {}
    for i in np.flatnonzero(effective):
        grouped.setdefault(source_targets[i], []).append(i)
    pairs = []
    for gene, indices in grouped.items():
        if len(indices) < 2:
            continue
        left, _ = response(indices[::2])
        right, _ = response(indices[1::2])
        mask = basal >= 5
        if gene in shared_positions:
            mask = mask.copy()
            mask[shared_positions[gene]] = False
        pairs.append((left, right, mask))
    records = []
    for alpha in (0., .15, .3, .5, .7, 1.):
        err, baseline, comparisons = 0., 0., 0
        for left, right, mask in pairs:
            err += float((np.square(alpha*left-right)+np.square(alpha*right-left))[mask].sum())
            baseline += float((left**2+right**2)[mask].sum())
            comparisons += 2*int(mask.sum())
        records.append({'shrinkage':alpha,'relative_lfc_mse':err/max(baseline,1e-12),
                        'n_gene_comparisons':comparisons,'n_construct_pairs':len(pairs)})
    selected = min(records,key=lambda r:r['relative_lfc_mse']) if pairs else None
    # Distinct promoter/transcript constructs are different interventions,
    # not technical replicates. Record their agreement diagnostic, but do not
    # treat it as cross-context validation or tune transfer from that proxy.
    transfer_alpha = .3
    effects = np.zeros((1,len(targets),len(shared)),dtype=np.float64)
    available = np.zeros((1,len(targets)),dtype=bool)
    counts = np.zeros((1,len(targets)),dtype=np.float64)
    coverage = []
    for t,target in enumerate(targets):
        indices = grouped.get(target, [])
        if indices:
            effects[0,t], counts[0,t] = response(indices)
            available[0,t] = True
        coverage.append({'target_gene':target,'n_reference_cell_lines':int(bool(indices)),
                         'K562_cells':int(counts[0,t]),'effective_source_constructs':len(indices),
                         'status':'transferred' if indices else 'control-based downstream fallback'})
    pd.DataFrame(records).to_csv(analysis/'reference_cross_validation.csv',index=False)
    pd.DataFrame(coverage).to_csv(analysis/'reference_target_coverage.csv',index=False)
    choice = {'shrinkage':transfer_alpha,'temperature':0.}
    summary = {'source_cell_lines':['K562'], 'source_genes':len(source_symbols),
               'unique_source_genes':len(set(source_symbols)), 'shared_genes':len(shared),
               'targets_with_reference':int(available.sum()),'targets_total':len(targets),
               'core_control_profiles':int(controls.sum()), 'source_control_cells':int(cells[controls].sum()),
               'source_construct_agreement_diagnostic':selected,'transfer_settings':choice,
               'transfer_shrinkage_assumption':'Fixed 0.3 log2-effect shrinkage after Poisson-noise soft thresholding; cross-context accuracy is unvalidated.',
               'validation_scope':'Exploratory source K562 construct agreement only; separate promoter interventions are not independent technical replicates. No cross-context validation or VCC score.',
               'candidates':records, 'source_file_sha256':sha256(source)}
    write_json(analysis/'reference_model.json',summary)
    np.savez_compressed(analysis/'reference_model_arrays.npz',effects=effects.astype(np.float32),
                        available=available,basal=basal[None,:].astype(np.float32),counts=counts,
                        genes=shared,targets=targets,cell_lines=np.array(['K562']),
                        dest_positions=dest_positions,shrinkage=transfer_alpha,temperature=0.)
    log(f"Reference coverage {summary['targets_with_reference']}/{len(targets)} targets, {len(shared)} genes; "
        f"source construct pairs {len(pairs)}; transfer shrinkage {transfer_alpha}")
    return effects, available, basal[None,:], dest_positions, choice, ['K562'], summary


def context_factors(mean_cpm, model, genes, targets, knockdown):
    effects, available, basal, positions, choice, lines, _ = model
    local = mean_cpm[positions]
    local = local*1e6/local.sum()
    correlations = np.array([np.corrcoef(np.log1p(local), np.log1p(b))[0, 1] for b in basal])
    w = np.exp(choice['temperature']*(correlations-1))
    weights = w[:, None]*available
    denominator = weights.sum(axis=0)
    lfc = np.einsum('lt,ltg->tg', weights, effects)/np.maximum(denominator[:, None], 1e-12)
    lfc *= choice['shrinkage']
    factors = np.ones((len(targets), len(genes)), dtype=np.float64)
    factors[:, positions] = np.exp2(lfc)
    gene_positions = {g: i for i, g in enumerate(genes)}
    for t, target in enumerate(targets):
        # The on-target response is biologically required but excluded from all
        # competition metrics. It never enters the source holdout loss above.
        factors[t, gene_positions[target]] = 1-knockdown
    diagnostics = {'source_profile_correlations': dict(zip(lines, correlations.tolist())),
                   'source_weights_before_per_target_renormalization': dict(zip(lines, (w/w.sum()).tolist())),
                   'downstream_lfc_quantiles': quantiles(lfc.ravel())}
    return factors, diagnostics


def simulate_counts(x, rows, factors, rng):
    """Thin decreased counts and add Poisson increments to increased counts."""
    pred = x[rows].copy()
    scale = factors[pred.indices]
    lower = scale < 1
    higher = scale > 1
    values = pred.data.astype(np.int64)
    values[lower] = rng.binomial(values[lower], scale[lower])
    values[higher] += rng.poisson(values[higher]*(scale[higher]-1))
    if values.max(initial=0) > np.iinfo(np.int32).max:
        raise ValueError('Individual predicted count overflows int32')
    pred.data = values.astype(np.int32)
    pred.eliminate_zeros()
    sums = np.asarray(pred.sum(axis=1, dtype=np.int64)).ravel()
    if (sums <= 0).any() or sums.max(initial=0) > MAX_COUNTS:
        raise ValueError('Predicted library size outside competition limits')
    return pred


class SparseH5Writer:
    """Append CSR blocks; use int64 global pointers, even beyond 2^31 entries."""
    def __init__(self, path, obs, genes):
        skeleton = ad.AnnData(X=sparse.csr_matrix((len(obs), len(genes)), dtype=np.int32),
                             obs=obs, var=pd.DataFrame(index=pd.Index(genes, dtype=object)))
        skeleton.write_h5ad(path)
        self.file = h5py.File(path, 'r+')
        del self.file['X']
        self.group = self.file.create_group('X')
        self.group.attrs.update({'encoding-type': 'csr_matrix', 'encoding-version': '0.1.0',
                                 'shape': np.array([len(obs), len(genes)], dtype=np.int64)})
        opts = dict(shape=(0,), maxshape=(None,), chunks=(1_048_576,), compression='gzip',
                    compression_opts=1, shuffle=True)
        self.data = self.group.create_dataset('data', dtype='int32', **opts)
        self.indices = self.group.create_dataset('indices', dtype='int32', **opts)
        self.indptr = self.group.create_dataset('indptr', shape=(len(obs)+1,), dtype='int64')
        self.indptr[0] = 0
        self.rows = 0
        self.nnz = 0
        self.expected_rows = len(obs)
    def append(self, block):
        block = block.tocsr()
        block.sum_duplicates()
        block.eliminate_zeros()
        block.sort_indices()
        end = self.nnz + block.nnz
        if end > MAX_NNZ or self.rows + block.shape[0] > self.expected_rows:
            raise ValueError('Prediction exceeds declared shape or density cap')
        self.data.resize((end,))
        self.indices.resize((end,))
        self.data[self.nnz:end] = block.data
        self.indices[self.nnz:end] = block.indices
        self.indptr[self.rows+1:self.rows+block.shape[0]+1] = block.indptr[1:].astype(np.int64)+self.nnz
        self.rows += block.shape[0]
        self.nnz = end
    def close(self):
        self.file.close()
        if self.rows != self.expected_rows:
            raise ValueError(f'Incomplete prediction: {self.rows}/{self.expected_rows} rows')


def validate_h5ad(path, genes, targets, contexts, n):
    obj = ad.read_h5ad(path, backed='r')
    assert obj.shape == (len(targets)*len(contexts)*n, len(genes))
    assert np.array_equal(obj.var_names.to_numpy(), genes)
    assert obj.obs_names.is_unique
    expected = {(c, t): n for c in contexts for t in targets}
    actual = obj.obs.groupby(['context', 'target_gene'], observed=True).size().to_dict()
    assert actual == expected
    obj.file.close()
    with h5py.File(path, 'r') as f:
        assert f['X'].attrs['encoding-type'] == 'csr_matrix'
        ptr = f['X/indptr'][:]
        assert ptr[0] == 0 and np.all(np.diff(ptr) > 0)
        assert int(ptr[-1]) == len(f['X/data']) == len(f['X/indices']) <= MAX_NNZ
        maximum = 0
        # Chunk on rows so every stored count and every library size is checked
        # without ever materializing the full prediction.
        for row in range(0, len(ptr)-1, 1000):
            stop = min(row+1000, len(ptr)-1)
            start_nnz, stop_nnz = int(ptr[row]), int(ptr[stop])
            data = f['X/data'][start_nnz:stop_nnz]
            indices = f['X/indices'][start_nnz:stop_nnz]
            assert np.issubdtype(data.dtype, np.integer) and np.all(data > 0)
            assert indices.min(initial=0) >= 0 and indices.max(initial=0) < len(genes)
            local_ptr = ptr[row:stop+1]-start_nnz
            counts = np.add.reduceat(data.astype(np.int64), local_ptr[:-1])
            maximum = max(maximum, int(counts.max()))
            assert counts.min() > 0 and counts.max() <= MAX_COUNTS
        return {'n_cells': len(ptr)-1, 'n_genes': len(genes), 'nnz': int(ptr[-1]),
                'maximum_cell_counts': maximum, 'groups': len(expected),
                'cells_per_group': n, 'all_stored_counts': 'positive integers',
                'gene_order': 'exact match', 'context_target_cell_counts': 'exact match'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--seed', type=int, default=20261003)
    parser.add_argument('--knockdown', type=float, default=.9)
    parser.add_argument('--analysis-only', action='store_true')
    parser.add_argument('--package-only', action='store_true')
    parser.add_argument('--force', action='store_true')
    args = parser.parse_args()
    root = args.root.resolve()
    if not 0 < args.knockdown < 1:
        parser.error('--knockdown must be strictly between 0 and 1')
    analysis = root/'analysis'
    analysis.mkdir(exist_ok=True)
    h5ad = root/'prediction.h5ad'
    package = root/'prediction.vcc'
    manifest, genes, targets = read_spec(root)
    contexts = manifest['contexts']
    n = manifest['cells_per_pert']
    start = time.time()
    if not args.package_only:
        if not args.analysis_only and not args.force and (h5ad.exists() or package.exists()):
            raise FileExistsError('Prediction already exists; use --force for an intentional rebuild')
        log('Inspecting and hashing the original project files')
        files = inventory(root, analysis)
        summaries, cpm_profiles, gene_frames = [], [], []
        model = None
        if not args.analysis_only:
            if not (root/'reference'/'K562_gwps_raw_bulk_01.h5ad').exists():
                subprocess.run([sys.executable, str(root/'download_reference.py'), '--root', str(root)], check=True)
            model = reference_model(root, genes, targets, analysis)
            context_labels = np.repeat(contexts, len(targets)*n)
            target_labels = np.tile(np.repeat(targets, n), len(contexts))
            obs = pd.DataFrame({'target_gene': pd.Categorical(target_labels, categories=targets),
                                'context': pd.Categorical(context_labels, categories=contexts)},
                               index=[f'pred_{i:06d}' for i in range(len(context_labels))])
            partial = h5ad.with_suffix('.h5ad.partial')
            writer = SparseH5Writer(partial, obs, genes)
        diagnostics, pred_records = {}, []
        for ci, context in enumerate(contexts):
            log(f'Analysing control context {context}')
            x, guides, libs, mean_cpm, summary, stats = load_control(root, context, genes, manifest)
            summaries.append(summary)
            cpm_profiles.append(mean_cpm)
            gene_frames.append(stats)
            stats.to_csv(analysis/f'gene_statistics_{context}.csv', index=False)
            if not args.analysis_only:
                factors, diag = context_factors(mean_cpm, model, genes, targets, args.knockdown)
                diagnostics[context] = diag
                for ti, target in enumerate(targets):
                    rng = np.random.default_rng(np.random.SeedSequence([args.seed, ci, ti]))
                    rows = guide_stratified_sample(guides, libs, n, rng)
                    block = simulate_counts(x, rows, factors[ti], rng)
                    writer.append(block)
                    target_index = int(np.flatnonzero(genes == target)[0])
                    source_mean = float(x[rows, target_index].mean())
                    predicted_mean = float(block[:, target_index].mean())
                    library_mean = float(block.sum(dtype=np.int64))/n
                    counts = Counter(guides[rows])
                    pred_records.append({'context': context, 'target_gene': target, 'n_cells': n,
                        'nnz': int(block.nnz), 'sampled_target_mean': source_mean,
                        'predicted_target_mean': predicted_mean,
                        'mean_predicted_library': library_mean,
                        'minimum_cells_per_control_guide': min(counts.values()),
                        'maximum_cells_per_control_guide': max(counts.values())})
                    if (ti+1) % 25 == 0:
                        log(f'Context {context}: {ti+1}/{len(targets)} targets; {writer.nnz:,} stored counts')
                del factors
            del x
            gc.collect()
        write_json(analysis/'control_summary.json', summaries)
        corr = np.corrcoef(np.log1p(cpm_profiles))
        pd.DataFrame(corr, index=contexts, columns=contexts).to_csv(analysis/'context_profile_correlations.csv')
        pd.concat(gene_frames).query('gene_name in @targets').to_csv(analysis/'target_gene_statistics.csv', index=False)
        write_json(analysis/'analysis_summary.json', {'original_files': len(files), 'manifest': manifest,
                    'context_profile_log_cpm_correlations': corr.tolist(),
                    'challenge_responses_available_locally': False})
        if args.analysis_only:
            log('Control analysis complete; no prediction requested')
            return
        writer.close()
        partial.replace(h5ad)
        pd.DataFrame(pred_records).to_csv(analysis/'prediction_group_statistics.csv', index=False)
        write_json(analysis/'context_transfer_diagnostics.json', diagnostics)
        versions = {name: importlib.metadata.version(name) for name in
                    ['vcc-cli', 'anndata', 'numpy', 'scipy', 'h5py', 'pandas', 'duckdb', 'pyarrow']}
        write_json(analysis/'run_parameters.json', {'created_utc': datetime.now(timezone.utc).isoformat(),
                    'seed': args.seed, 'knockdown_fraction_assumption': args.knockdown,
                    'model': 'Genome-wide K562 CRISPR response transfer with count-noise suppression and fixed 0.3 shrinkage',
                    'versions': versions, 'script_sha256': sha256(Path(__file__)),
                    'limitations': ['No challenge perturbation ground truth: no local VCC score.',
                                   'Transferred effects cover only reference-measured genes and targets.',
                                   'No induction in genes absent from all context controls.',
                                   'Cell-line and assay differences can limit transfer.']})
    log('Verifying every stored count, every cell library, labels, shape and gene order')
    verified = validate_h5ad(h5ad, genes, targets, contexts, n)
    write_json(analysis/'prediction_validation.json', verified)
    log('Running official vcc prep validation and packaging')
    partial_package = package.with_suffix('.vcc.partial')
    command = [sys.executable, '-m', 'vcc', '--json', 'prep', str(h5ad),
               '-g', str(root/'gene_names.csv'), '--perts', str(root/'pert_counts.csv'),
               '--contexts', ','.join(contexts), '--cells-per-pert', str(n),
               '-o', str(partial_package)]
    if args.force:
        command += ['--force']
    result = subprocess.run(command, capture_output=True, text=True, encoding='utf-8', errors='replace')
    (analysis/'official_prep_stdout.json').write_text(result.stdout, encoding='utf-8')
    (analysis/'official_prep_stderr.txt').write_text(result.stderr, encoding='utf-8')
    if result.returncode:
        raise RuntimeError(f'Official prep failed; see analysis/official_prep_stderr.txt\n{result.stderr[-3000:]}')
    from vcc.vccfile import validate_vcc, read_vcc_meta
    validate_vcc(partial_package)
    meta = read_vcc_meta(partial_package)
    assert meta and meta['nnz'] == verified['nnz']
    partial_package.replace(package)
    write_json(analysis/'submission_ready.json', {'file': str(package), 'bytes': package.stat().st_size,
               'sha256': sha256(package), 'official_prep_passed': True, 'package_meta': meta,
               'prediction_validation': verified, 'elapsed_seconds': time.time()-start,
               'leaderboard_score': None, 'uploaded': False})
    log(f'READY: {package} ({package.stat().st_size/1e9:.3f} GB); not uploaded')


if __name__ == '__main__':
    main()
