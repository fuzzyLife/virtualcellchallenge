"""Independent verification of identifiers, labels, metrics and raw-rank spot checks."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from threadpoolctl import threadpool_limits
from analyze_target_rank_correlations import ROOT, digest

OUT = ROOT / 'analysis/full_rank_clustering'


def run():
    report = json.loads((OUT / 'analysis_report.json').read_text())
    for name, record in report['outputs'].items():
        assert digest(OUT / name) == record['sha256'], name
    for name, sha in report['scripts'].items():
        assert digest(ROOT / name) == sha, name
    cells = pd.read_csv(OUT / 'per_cell_full_rank_clusters.csv')
    genes = pd.read_csv(OUT / 'gene_axis.csv').gene_name.to_numpy()
    np.testing.assert_array_equal(genes, pd.read_csv(ROOT / 'gene_names.csv').gene_name)
    assert cells.shape[0] == 55200 and len(genes) == 18533
    assert not cells.cell_id.duplicated().any()
    assert not cells.select_dtypes(include='number').isna().any().any()
    assert not report['labels_used_for_fit'] and not report['pca_used_for_fit']
    models = np.load(OUT / 'cluster_models.npz', allow_pickle=False)
    metrics = pd.read_csv(OUT / 'cluster_metrics.csv')
    for row in metrics.itertuples():
        key = f'k{row.k}_repeat{row.repeat}'
        label = cells['cluster_' + (f'k{row.k}' if row.repeat == 0 else key)]
        assert set(label) == set(range(row.k))
        np.testing.assert_allclose(np.linalg.norm(models[key], axis=1), 1, atol=1e-6)
        counts = pd.crosstab(cells.context, label)
        assert abs(row.context_purity - counts.max(axis=0).sum() / len(cells)) < 1e-12
        assert abs(row.context_ARI - adjusted_rand_score(cells.context, label)) < 1e-12
        assert abs(row.context_NMI - normalized_mutual_info_score(cells.context, label)) < 1e-12
        assert row.min_cluster_cells == counts.sum(axis=0).min()
        assert row.singleton_clusters == (counts.sum(axis=0) == 1).sum()
        sample = pd.read_csv(OUT / 'silhouette_sample.csv').sample_index.to_numpy()
        assert row.sample_clusters_represented == len(np.unique(label.iloc[sample]))
        assert row.converged
    composition = pd.read_csv(OUT / 'cluster_composition_k3.csv')
    expected = pd.crosstab(cells.context, cells.cluster_k3)
    actual = composition.pivot(index='context', columns='cluster', values='cells')
    np.testing.assert_array_equal(actual, expected)
    assert abs(report['k3_repeat_ARI'] - adjusted_rand_score(cells.cluster_k3, cells.cluster_k3_repeat1)) < 1e-12
    for name, groupcols in [('depth_quintile_composition.csv', ['context']),
                            ('guide_composition.csv', ['context'])]:
        frame = pd.read_csv(OUT / name)
        assert frame.cells.sum() == 55200
        assert (frame.groupby(groupcols).cells.sum() == 18400).all()
    rng = np.random.default_rng(93)
    raw_checks = []
    sample = pd.read_csv(OUT / 'silhouette_sample.csv')
    sample_ranks = np.empty((len(sample), len(genes)), dtype=np.float64)
    for context in ['A', 'B', 'C']:
        record = next(x for x in report['inputs'] if x['file'] == f'context_{context}.h5ad')
        assert digest(ROOT / record['file']) == record['sha256']
        obj = ad.read_h5ad(ROOT / record['file'], backed='r')
        try:
            sub = cells[cells.context == context].reset_index(drop=True)
            np.testing.assert_array_equal(sub.cell_id, obj.obs_names)
            np.testing.assert_array_equal(sub.ntc_id, obj.obs.ntc_id.astype(str))
            np.testing.assert_array_equal(obj.var_names, genes)
            indices = np.sort(rng.choice(len(obj), 20, replace=False))
            raw = obj.X[indices].toarray()
            # Independently rank each original row in double precision.
            z = np.array([rankdata(row, method='average') for row in raw])
            z -= z.mean(axis=1, keepdims=True)
            z /= np.linalg.norm(z, axis=1, keepdims=True)
            np.testing.assert_array_equal(sub.total_counts.iloc[indices], raw.sum(axis=1))
            np.testing.assert_array_equal(sub.detected_genes.iloc[indices], (raw > 0).sum(axis=1))
            for row in metrics.itertuples():
                key = f'k{row.k}_repeat{row.repeat}'
                column = 'cluster_' + (f'k{row.k}' if row.repeat == 0 else key)
                np.testing.assert_array_equal((z @ models[key].T).argmax(axis=1), sub[column].iloc[indices])
            similarity = z @ models['k3_repeat0'].T
            np.testing.assert_allclose(similarity.max(axis=1), sub.own_cluster_centroid_rho.iloc[indices], atol=2e-6)
            np.testing.assert_allclose((z - models['pca_mean']) @ models['pca_components'].T,
                                       sub.iloc[indices][[f'PC{i}' for i in range(1, 9)]], atol=2e-5)
            raw_checks.extend(raw[:5])
            membership = sample.cell_id.isin(sub.cell_id).to_numpy()
            positions = sample.loc[membership, 'sample_index'].to_numpy() - ['A','B','C'].index(context)*18400
            sample_raw = obj.X[positions].toarray()
            ranks = rankdata(sample_raw, axis=1, method='average')
            ranks -= ranks.mean(axis=1, keepdims=True)
            ranks /= np.linalg.norm(ranks, axis=1, keepdims=True)
            sample_ranks[membership] = ranks
        finally:
            obj.file.close()
    for i in range(len(raw_checks)):
        for j in range(i):
            a = rankdata(raw_checks[i]); b = rankdata(raw_checks[j])
            a -= a.mean(); b -= b.mean()
            value = (a @ b) / (np.linalg.norm(a) * np.linalg.norm(b))
            assert abs(value - spearmanr(raw_checks[i], raw_checks[j]).statistic) < 1e-12
    with threadpool_limits(limits=12):
        distances = np.clip(1 - sample_ranks @ sample_ranks.T, 0, 2)
        np.fill_diagonal(distances, 0)
    for row in metrics.itertuples():
        column = 'cluster_' + (f'k{row.k}' if row.repeat == 0 else f'k{row.k}_repeat{row.repeat}')
        value = silhouette_score(distances, cells[column].iloc[sample.sample_index], metric='precomputed')
        assert abs(value - row.sample_silhouette_spearman) < 2e-6
    summary = {'status': 'passed', 'all_cell_identifier_checks': 55200,
               'original_full_gene_cells_spot_checked': 60, 'cluster_assignments_spot_checked': 360,
               'scipy_spearman_cell_pair_checks': 105,
               'silhouette_sample_independently_reranked_cells': 1500, 'silhouette_scores_checked': 6,
               'checks': 'Input/output/script SHA256, gene axis, all IDs, compositions, ARI/NMI/purity, all fits converged, six-model assignments and PCA from raw rows'}
    (OUT / 'validation.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    run()
