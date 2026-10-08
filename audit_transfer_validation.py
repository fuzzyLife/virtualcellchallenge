"""Audit frozen V7 validation against common-control VCC2026 profile semantics.

This is a retrospective diagnostic, not a new model-selection benchmark. It
verifies saved profile hashes and reconstructs the disjoint control populations.
No V7 numerical files, selections, or submission payloads are changed.
"""
from pathlib import Path
import argparse
import gc
import hashlib
import json
import platform
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parent
OLD = ROOT / 'experiments/v7_multisource_adaptive/analysis'
OUT = ROOT / 'analysis/transfer_validation_audit'
CANDIDATES = ['v4_transfer', 'coverage_fill', 'coverage_uniform', 'calibrated_pool',
              'pool_low_noise', 'pool_cal10', 'neural_response', 'neural_blend25', 'ridge_blend25']
BASE_SEED = 20261008


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def pds(prediction, reference, control):
    """Raw PDS in a fixed feature space with exactly one context control."""
    prediction, reference, control = map(lambda a: np.asarray(a, dtype=np.float64),
                                         (prediction, reference, control))
    if prediction.ndim != 2 or reference.shape != prediction.shape:
        raise ValueError('Prediction and reference must have the same 2D shape')
    if control.shape != (prediction.shape[1],):
        raise ValueError('Control must be one common 1D profile, not a profile per target')
    if len(prediction) < 2 or not prediction.shape[1]:
        raise ValueError('At least two targets and one scored gene are required')
    if not all(np.isfinite(a).all() for a in (prediction, reference, control)):
        raise ValueError('Profiles must be finite')
    return distance_rank(prediction-control, reference-control)


def distance_rank(a, b):
    norm = np.linalg.norm(a, axis=1)[:, None] * np.linalg.norm(b, axis=1)[None, :]
    distance = 1 - np.divide(a @ b.T, norm, out=np.zeros((len(a), len(a))), where=norm > 0)
    diagonal = np.diag(distance)
    rank = (distance < diagonal[:, None]).sum(axis=1)
    rank = rank + .5 * ((distance == diagonal[:, None]).sum(axis=1)-1)
    return float((1-rank/(len(a)-1)).mean())


def corrected_expression(pred, truth, ctrl, pj, tj, cj, own):
    """Documented profile ratio, whole-gene scalar jackknife corrections.

Retains signed numerator. Reports an unresolved denominator rather than silently
replacing it with 1e-12. No leaderboard anchors are estimated here.
"""
    mask = np.ones(pred.shape, dtype=bool)
    mask[np.arange(len(own)), own] = False
    sizes = mask.sum(axis=1)
    pc, tc = pj.sum(axis=1)/sizes, tj.sum(axis=1)/sizes
    cc = float(cj.sum())/sizes
    weights = mask/sizes[:, None]
    mean = np.divide((pred*weights).sum(axis=0), weights.sum(axis=0),
                     out=np.zeros(pred.shape[1]), where=weights.sum(axis=0) > 0)
    spread = float(((pred-mean)**2*weights).sum())
    cap = np.minimum(pc, tc)
    rho = min(1., spread/cap.sum()) if cap.sum() > 0 else 1.
    numerator = float((((pred-truth)**2*mask).sum(axis=1)/sizes-rho*cap-tc).sum())
    denominator = float((((truth-ctrl)**2*mask).sum(axis=1)/sizes-tc-cc).sum())
    return {'common_expression_numerator': numerator,
            'common_expression_denominator': denominator,
            'common_expression_ratio': numerator/denominator if denominator > 0 else None,
            'common_expression_denominator_resolved': bool(denominator > 0),
            'common_expression_rho': float(rho)}


def profile(x):
    summed = np.asarray(sparse.csr_matrix(x, dtype=np.float64).sum(axis=0)).ravel()
    if summed.sum() <= 0:
        raise ValueError('Empty control population')
    return np.log1p(50000*summed/summed.sum())


def profile_jackknife_chunked(x, chunk=256):
    """Delete-one jackknife using mergeable central moments (bounded dense RAM)."""
    x = sparse.csr_matrix(x, dtype=np.float64)
    counts = np.asarray(x.sum(axis=0)).ravel()
    depth = np.asarray(x.sum(axis=1)).ravel()
    total = depth.sum()
    result = np.log1p(50000*counts/total)
    if len(depth) < 2:
        return result, np.zeros_like(result)
    if (depth <= 0).any() or (total-depth <= 0).any():
        raise ValueError('Invalid library in jackknife')
    n = 0
    mean = np.zeros_like(result)
    m2 = np.zeros_like(result)
    for start in range(0, len(depth), chunk):
        block = x[start:start+chunk].toarray()
        deleted = np.log1p(50000*(counts-block)/(total-depth[start:start+chunk])[:, None])
        bn = len(block)
        bm = deleted.mean(axis=0)
        delta = bm-mean
        m2 += ((deleted-bm)**2).sum(axis=0) + delta**2*(n*bn/(n+bn))
        mean += delta*(bn/(n+bn))
        n += bn
    return result, (n-1)*m2/n


def self_test():
    rng = np.random.default_rng(31)
    truth = rng.normal(size=(9, 21))
    ctrl = rng.normal(size=21)
    constant = np.broadcast_to(rng.normal(size=21), truth.shape)
    np.testing.assert_allclose(pds(constant, truth, ctrl), .5, atol=1e-14)
    np.testing.assert_allclose(pds(np.broadcast_to(ctrl, truth.shape), truth, ctrl), .5)
    np.testing.assert_allclose(pds(truth, truth, ctrl), 1.)
    perm = rng.permutation(len(truth))
    np.testing.assert_allclose(pds(constant[perm], truth[perm], ctrl), .5)
    try:
        pds(truth, truth, np.broadcast_to(ctrl, truth.shape))
    except ValueError:
        pass
    else:
        raise AssertionError('Target-specific controls were not rejected')
    # Historical semantics give a strong score with no target-specific prediction.
    variable_controls = rng.normal(size=truth.shape)*20
    bad = distance_rank(constant-variable_controls, truth-variable_controls)
    assert bad > .9, bad
    # Zero-error and control-pasting limits with zero sampling corrections.
    zeros = np.zeros_like(truth)
    own = np.arange(len(truth))
    perfect = corrected_expression(truth, truth, ctrl, zeros, zeros, np.zeros_like(ctrl), own)
    null = corrected_expression(np.broadcast_to(ctrl, truth.shape), truth, ctrl, zeros, zeros, np.zeros_like(ctrl), own)
    np.testing.assert_allclose(perfect['common_expression_ratio'], 0.)
    np.testing.assert_allclose(null['common_expression_ratio'], 1.)
    unresolved = corrected_expression(truth, truth, ctrl, zeros, np.full_like(truth, 100.), np.zeros_like(ctrl), own)
    assert unresolved['common_expression_ratio'] is None
    assert not unresolved['common_expression_denominator_resolved']
    from v6_methods import profile_jackknife
    x = sparse.csr_matrix(rng.poisson(4, size=(17, 21)))
    expected = profile_jackknife(x)
    actual = profile_jackknife_chunked(x, chunk=5)
    for a, b in zip(actual, expected):
        np.testing.assert_allclose(a, b, atol=1e-12, rtol=1e-10)
    return {'constant_prediction_pds': .5, 'zero_effect_pds': .5,
            'perfect_prediction_pds': 1., 'historical_synthetic_constant_pds': bad,
            'common_control_shape_guard': 'passed', 'chunked_jackknife_equivalence': 'passed',
            'expression_limits_and_denominator_guard': 'passed'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    tests = self_test()
    if args.self_test:
        print(json.dumps(tests, indent=2))
        return
    # Imported only for full audit; this loader verifies all extracted count hashes.
    from fit_v7_response import raw_inputs
    OUT.mkdir(parents=True, exist_ok=True)
    axis_hash = sha256(OLD/'source_cache.npz')
    axis_manifest = json.loads((OLD/'source_provenance.json').read_text())
    if axis_hash != axis_manifest['cache_sha256']:
        raise ValueError('Source axis cache hash mismatch')
    with np.load(OLD/'source_cache.npz', allow_pickle=False) as z:
        genes, targets = z['genes'], z['targets']
    gene_lookup = {str(g): i for i, g in enumerate(genes)}
    scored = ~np.isin(genes, targets)
    rows, provenance, populations = [], {}, {}
    for family in ['H1', 'HCT116', 'HEK293T']:
        print('Reading verified source', family, flush=True)
        x, meta, measured, source = raw_inputs(family, genes)
        controls = np.flatnonzero(meta.gene_target.to_numpy() == 'Non-Targeting')
        rng = np.random.default_rng(BASE_SEED)
        left, right = [], []
        for _, group in meta.iloc[controls].groupby('sample', sort=True):
            ids = rng.permutation(group.index.to_numpy())
            half = len(ids)//2
            left.extend(ids[:half]); right.extend(ids[half:])
        assert not set(left) & set(right)
        ctrl, cj = profile_jackknife_chunked(x[right])
        numerical_dependencies = {f: sha256(ROOT/f) for f in ['fit_v7_response.py', 'v6_methods.py', 'build_submission.py']}
        constant = profile(x[left])
        id_lookup = {str(cid): i for i, cid in enumerate(meta.cell_id)}
        left_ids, right_ids = set(meta.iloc[left].cell_id.astype(str)), set(meta.iloc[right].cell_id.astype(str))
        records = []
        for file in sorted((OLD/'cell_evaluation').glob(f'{family}_*.npz')):
            side = file.with_suffix('.json')
            rec = json.loads(side.read_text())
            digest = sha256(file)
            if digest != rec['sha256']:
                raise ValueError(f'Frozen profile hash mismatch: {file.name}')
            cells = rec['cells']
            template_ids = list(map(str, cells['template']))
            reference_ids = list(map(str, cells['reference']))
            assert set(template_ids) <= left_ids and set(reference_ids) <= right_ids
            assert not set(template_ids) & set(reference_ids)
            with np.load(file, allow_pickle=False) as z:
                arrays = {k: z[k].copy() for k in z.files}
            arrays['template_only'], arrays['template_only_jack'] = profile_jackknife_chunked(
                x[[id_lookup[c] for c in template_ids]])
            records.append((cells, arrays))
            provenance[file.name] = {'sha256': digest, 'sidecar_sha256': sha256(side),
                                     'frozen_inputs': rec['inputs']}
        for seed in sorted({r[0]['seed'] for r in records}):
            for split in ['calibration', 'confirmation', 'all']:
                selected = [(c, a) for c, a in records if c['seed'] == seed and
                            (split == 'all' or c['split'] == split)]
                names = [c['target_gene'] for c, _ in selected]
                truth = np.array([a['truth'] for _, a in selected])
                matched = np.array([a['ctrl'] for _, a in selected])
                tj = np.array([a['truth_jack'] for _, a in selected])
                own = np.array([gene_lookup[t] for t in names])
                for candidate in CANDIDATES + ['constant_control', 'template_only']:
                    if candidate == 'constant_control':
                        pred = np.broadcast_to(constant, truth.shape)
                    elif candidate == 'template_only':
                        pred = np.array([a[candidate] for _, a in selected])
                    else:
                        pred = np.array([a[candidate+'_pred'] for _, a in selected])
                    old_pds = distance_rank((pred-matched)[:, scored], (truth-matched)[:, scored])
                    common_pds = pds(pred[:, scored], truth[:, scored], ctrl[scored])
                    if candidate == 'constant_control':
                        np.testing.assert_allclose(common_pds, .5, atol=1e-12)
                    row = {'family': family, 'seed': seed, 'split': split, 'candidate': candidate,
                           'targets': len(names), 'historical_pds': old_pds, 'common_control_pds': common_pds,
                           'pds_change': common_pds-old_pds}
                    if candidate != 'constant_control':
                        pj = np.array([a[candidate+'_jack'] for _, a in selected])
                        row.update(corrected_expression(pred, truth, ctrl, pj, tj, cj, own))
                    rows.append(row)
        populations[family] = {'template_population_cells': len(left), 'reference_population_cells': len(right),
                                'control_batches': int(meta.iloc[controls]['sample'].nunique()),
                                'numerical_dependencies_sha256': numerical_dependencies,
                                'reference_cell_ids': meta.iloc[right].cell_id.tolist(),
                                'template_cell_ids': meta.iloc[left].cell_id.tolist(), 'source': source}
        pd.DataFrame(rows).to_csv(OUT/'profile_audit.csv', index=False)
        del x, meta, records
        gc.collect()
        print('Audited', family, 'with one shared reference control', flush=True)
    frame = pd.DataFrame(rows)
    historical = pd.concat([pd.read_csv(OLD/f'cell_summary_{line}.csv') for line in ['H1', 'HCT116', 'HEK293T']])
    comparison = frame.merge(historical[['family','seed','split','candidate','pds']],
                             on=['family','seed','split','candidate'], how='inner', validate='one_to_one')
    assert len(comparison) == 162
    np.testing.assert_allclose(comparison.historical_pds, comparison.pds, atol=1e-12, rtol=0)
    summary = frame.groupby(['family','split','candidate'], sort=True).agg(
        seeds=('seed','size'), targets=('targets','first'), historical_pds=('historical_pds','mean'),
        common_control_pds=('common_control_pds','mean'),
        common_expression_ratio=('common_expression_ratio','mean')).reset_index()
    summary.to_csv(OUT/'profile_summary.csv', index=False)
    nulls = frame[frame.candidate == 'constant_control']
    result = {'purpose': 'Retrospective metric-origin and null audit; not model selection or leaderboard estimates',
              'self_tests': tests, 'verified_profile_files': len(provenance),
              'historical_pds_reproduced_rows': len(comparison), 'source_axis_sha256': axis_hash,
              'constant_historical_pds_range': [float(nulls.historical_pds.min()), float(nulls.historical_pds.max())],
              'constant_common_pds_range': [float(nulls.common_control_pds.min()), float(nulls.common_control_pds.max())],
              'evaluation_caveats': [
                  'Held public panels contain 11/18/18 targets, not the official 300-target contexts or anchors.',
                  'Common source controls pool different technical batches; this is not a replacement biological validation panel.',
                  'Saved templates are matched to actual treated batch proportions, information unavailable for challenge targets.',
                  'Historical confirmation data were already inspected. No tuning, reselection or submission modification occurs.',
                  'DE metrics are not recomputed: saved profiles cannot reconstruct per-cell rank tests.'],
              'script_sha256': sha256(__file__), 'python': platform.python_version(), 'numpy': np.__version__,
              'metric_source': 'https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics.md',
              'frozen_profile_provenance': provenance, 'control_populations': populations,
              'output_sha256': sha256(OUT/'profile_audit.csv'),
              'summary_sha256': sha256(OUT/'profile_summary.csv')}
    (OUT/'audit_manifest.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(frame[(frame.split == 'all') & frame.candidate.isin(['constant_control', 'template_only', 'v4_transfer'])]
          [['family','seed','candidate','historical_pds','common_control_pds']].to_string(index=False))
    print('Saved', OUT, flush=True)


if __name__ == '__main__':
    with threadpool_limits(limits=4):
        main()
