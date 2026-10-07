"""Compare frozen V4/V5 effect changes to inform a V6 plan, without fitting or uploading."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from build_submission import ROOT, sha256, write_json

OUT = ROOT / 'analysis/v6_strategy'


def cosine(a, b):
    norm = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / norm) if norm > 0 else None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    paths = {v: ROOT / f'experiments/{folder}/analysis/reference_model_arrays.npz'
        for v, folder in [('v4', 'v4_lineage_response'), ('v5', 'v5_atlas_calibrated')]}
    with np.load(paths['v4'], allow_pickle=False) as z:
        old, genes, targets = z['effects'].copy(), z['genes'].copy(), z['targets'].copy()
    with np.load(paths['v5'], allow_pickle=False) as z:
        new = z['effects'].copy()
        np.testing.assert_array_equal(z['genes'], genes)
        np.testing.assert_array_equal(z['targets'], targets)
    assert old.shape == new.shape == (3, 300, 18533)
    assert np.isfinite(old).all() and np.isfinite(new).all()
    panel = np.isin(genes, targets)
    rows = []
    for ci, context in enumerate('ABC'):
        stats = pd.read_csv(ROOT / f'analysis/gene_statistics_{context}.csv')
        np.testing.assert_array_equal(stats.gene_name, genes)
        expressed = stats.mean_cpm.to_numpy() > 5
        for ti, target in enumerate(targets):
            mask = expressed & (genes != target)
            a, b = old[ci, ti, mask].astype(float), new[ci, ti, mask].astype(float)
            both = (np.abs(a) >= .05) & (np.abs(b) >= .05)
            delta = new[ci, ti].astype(float) - old[ci, ti].astype(float)
            delta[genes == target] = 0
            rows.append({'context': context, 'target_gene': target,
                'changed_output_coordinates': int(np.count_nonzero(np.abs(delta) > 1e-7)),
                'expressed_changed_coordinates': int(np.count_nonzero(np.abs(b-a) > 1e-7)),
                'both_effect_magnitude_at_least_0_05': int(both.sum()),
                'opposite_signs_among_both': int(np.count_nonzero(a[both] * b[both] < 0)),
                'v4_active_expressed_genes': int(np.count_nonzero(np.abs(a) >= .05)),
                'v5_active_expressed_genes': int(np.count_nonzero(np.abs(b) >= .05)),
                'mean_absolute_expressed_effect_change': float(np.mean(np.abs(b-a))),
                'expressed_effect_cosine': cosine(a,b),
                'nonpanel_effect_cosine': cosine(old[ci,ti,~panel].astype(float), new[ci,ti,~panel].astype(float))})
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / 'v4_v5_target_effect_changes.csv', index=False)
    summary = []
    for context, group in frame.groupby('context'):
        count = int(group.both_effect_magnitude_at_least_0_05.sum())
        summary.append({'context': context, 'changed_target_groups': int((group.changed_output_coordinates > 0).sum()),
            'changed_output_coordinates': int(group.changed_output_coordinates.sum()),
            'opposite_signs_both_active': int(group.opposite_signs_among_both.sum()),
            'opposite_sign_fraction_both_active': float(group.opposite_signs_among_both.sum()/count) if count else None,
            'mean_v4_active_expressed_genes': float(group.v4_active_expressed_genes.mean()),
            'mean_v5_active_expressed_genes': float(group.v5_active_expressed_genes.mean()),
            'median_expressed_effect_cosine': float(group.expressed_effect_cosine.median()),
            'median_nonpanel_effect_cosine': float(group.nonpanel_effect_cosine.median())})
    pd.DataFrame(summary).to_csv(OUT / 'v4_v5_context_effect_changes.csv', index=False)
    write_json(OUT / 'effect_comparison_scope.json', {'source_sha256': {k:sha256(v) for k,v in paths.items()},
        'script_sha256': sha256(Path(__file__)), 'contexts': summary, 'groups': 900,
        'scope': 'Model mean-response log2FC arrays only; not generated counts or measured perturbation truth.',
        'eligible_genes': 'Provided context control mean CPM >5, excluding own target; nonpanel cosine excludes all 300 inputs.',
        'effect_threshold': '0.05 absolute log2FC is descriptive, not a significant-call threshold.',
        'limitation': 'Sign changes between models do not identify incorrect signs or explain the published fidelity loss.'})
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == '__main__':
    main()
