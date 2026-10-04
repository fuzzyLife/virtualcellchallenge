"""Unsupervised spherical clustering of all VCC controls in full Spearman geometry."""
from pathlib import Path
from datetime import datetime, timezone
from importlib.metadata import version
import argparse
import json
import platform
import time
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.stats import rankdata, spearmanr
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_target_rank_correlations import ROOT, digest, table

OUT = ROOT / 'analysis/full_rank_clustering'
CONTEXTS = ['A', 'B', 'C']


def rank_unit_rows(raw):
    """Exact complete-data Spearman features; reject gaps rather than impute them."""
    raw = np.asarray(raw)
    if raw.ndim != 2 or not np.isfinite(raw).all():
        raise ValueError('A common, entirely finite gene axis is required')
    z = rankdata(raw, method='average', axis=1)
    z -= z.mean(axis=1, keepdims=True)
    norm = np.linalg.norm(z, axis=1)
    if np.any(norm == 0):
        raise ValueError('Constant cell vectors have undefined rank correlations')
    return np.asarray(z / norm[:, None], dtype=np.float32, order='C')


def spherical_refine(z, centers, max_iter=100):
    """Maximize sum of cell-to-unit-centroid correlations, on every full gene row."""
    centers = np.array(centers, dtype=np.float32, copy=True)
    norms = np.linalg.norm(centers, axis=1)
    if np.any(norms == 0):
        raise ValueError('Undefined initial centroid')
    centers /= norms[:, None]
    previous = None
    history = []
    for iteration in range(max_iter):
        similarity = z @ centers.T
        labels = similarity.argmax(axis=1)
        counts = np.bincount(labels, minlength=len(centers))
        if np.any(counts == 0):
            raise ValueError('Empty spherical cluster')
        changed = len(z) if previous is None else int(np.sum(labels != previous))
        history.append({'iteration': iteration, 'changed_cells': changed,
                        'sum_centroid_rho': float(similarity[np.arange(len(z)), labels].sum(dtype=np.float64))})
        if changed == 0:
            return labels, centers, history, True
        membership = sparse.csr_matrix((np.ones(len(z), np.float32), (labels, np.arange(len(z)))),
                                      shape=(len(centers), len(z)))
        centers = np.asarray(membership @ z)
        centers /= np.linalg.norm(centers, axis=1, keepdims=True)
        previous = labels
    # Final assignments always correspond to the saved centers.
    labels = (z @ centers.T).argmax(axis=1)
    return labels, centers, history, False


def cluster(z, k, seed):
    initial = MiniBatchKMeans(n_clusters=k, random_state=seed, n_init=3,
                             batch_size=1024, init_size=8192, max_iter=30,
                             max_no_improvement=15, reassignment_ratio=0)
    initial.fit(z)
    return spherical_refine(z, initial.cluster_centers_)


def separation(z, labels, distance):
    return float(silhouette_score(distance, labels, metric='precomputed'))


def make_plot(cells, composition, metrics, explained):
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), layout='constrained')
    order = np.random.default_rng(20261004).permutation(len(cells))
    x = cells.iloc[order]
    colors = ['#4477aa', '#ee7733', '#228833']
    for context, color in zip(CONTEXTS, colors):
        sub = x[x.context == context]
        axes[0, 0].scatter(sub.PC1, sub.PC2, s=1, alpha=.25, c=color, rasterized=True, label=context)
    axes[0, 0].legend(markerscale=7, title='Provided context')
    axes[0, 0].set_title('All cells, colored by A/B/C')
    for k, color in enumerate(colors):
        sub = x[x.cluster_k3 == k]
        axes[0, 1].scatter(sub.PC1, sub.PC2, c=color, s=1, alpha=.3, rasterized=True, label=str(k))
    axes[0, 1].legend(markerscale=7, title='Blind cluster')
    axes[0, 1].set_title('Blind full-matrix clustering: K=3')
    image = axes[0, 2].scatter(x.PC1, x.PC2, c=np.log10(x.total_counts), cmap='viridis', s=1, alpha=.3, rasterized=True)
    fig.colorbar(image, ax=axes[0, 2], label='log10(total raw counts)')
    axes[0, 2].set_title('Same projection, colored by count depth')
    for ax in axes[0]:
        ax.set_xlabel(f'PC1 ({explained[0]:.1%} of rank variance)')
        ax.set_ylabel(f'PC2 ({explained[1]:.1%} of rank variance)')
    matrix = composition.pivot(index='context', columns='cluster', values='cells').reindex(CONTEXTS)
    axes[1, 0].imshow(matrix, cmap='Blues', aspect='auto')
    for i in range(3):
        for j in range(3):
            axes[1, 0].text(j, i, f'{matrix.iloc[i,j]:,}', ha='center', va='center',
                            color='white' if matrix.iloc[i,j] > matrix.to_numpy().max()/2 else 'black')
    axes[1, 0].set_xticks(range(3), matrix.columns)
    axes[1, 0].set_yticks(range(3), CONTEXTS)
    axes[1, 0].set_xlabel('Cluster number (arbitrary)'); axes[1, 0].set_ylabel('Context')
    axes[1, 0].set_title('All 55,200 cells: context versus cluster')
    main = metrics[metrics.repeat == 0]
    axes[1, 1].plot(main.k, main.sample_silhouette_spearman, 'o-', label='Blind clusters')
    axes[1, 1].set_xticks(main.k); axes[1, 1].set_xlabel('Number of clusters')
    axes[1, 1].set_ylabel('Silhouette on 1 - Spearman rho')
    axes[1, 1].set_title('Exact full-gene distances; fixed 1,500-cell sample')
    axes[1, 2].boxplot([cells.loc[cells.cluster_k3 == k, 'total_counts'] for k in range(3)],
                       tick_labels=list(range(3)), showfliers=False)
    axes[1, 2].set_xlabel('Cluster'); axes[1, 2].set_ylabel('Total raw counts')
    axes[1, 2].set_title('Count-depth dispersion within clusters')
    fig.suptitle('Complete VCC control matrix: 55,200 cells x 18,533 genes\nRaw Spearman ranks; zeros retained; labels used only after fitting', fontsize=14)
    fig.savefig(OUT / 'full_rank_clusters.png', dpi=180)
    fig.savefig(OUT / 'full_rank_clusters.pdf')
    plt.close(fig)


def run(seed=20261004, threads=12):
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(ROOT / 'gene_names.csv').gene_name.to_numpy()
    if len(genes) != 18533 or len(set(genes)) != len(genes):
        raise ValueError('Unexpected gene axis')
    z = np.empty((55200, len(genes)), dtype=np.float32)
    obs = []
    inputs = []
    offset = 0
    old = json.loads((ROOT / 'analysis/target_rank_correlations/analysis_report.json').read_text())
    for context in CONTEXTS:
        path = ROOT / f'context_{context}.h5ad'
        sha = digest(path)
        expected = next(x for x in old['competition_inputs'] if x['file'] == path.name)
        if sha != expected['sha256']:
            raise ValueError('Original competition data changed')
        obj = ad.read_h5ad(path, backed='r')
        try:
            if obj.shape != (18400, len(genes)) or not np.array_equal(obj.var_names, genes):
                raise ValueError('Input shape or gene order differs')
            if set(obj.obs.context.astype(str)) != {context} or set(obj.obs.target_gene.astype(str)) != {'non-targeting'}:
                raise ValueError('Unexpected context or perturbation')
            frame = obj.obs.copy()
            frame.insert(0, 'cell_id', obj.obs_names)
            totals = []; detected = []
            for start in range(0, len(obj), 512):
                end = min(start + 512, len(obj))
                raw = obj.X[start:end].toarray()
                if np.any(raw < 0) or not np.equal(raw, np.floor(raw)).all():
                    raise ValueError('Expected raw nonnegative integer counts')
                z[offset + start:offset + end] = rank_unit_rows(raw)
                totals.extend(raw.sum(axis=1, dtype=np.float64))
                detected.extend((raw > 0).sum(axis=1))
            frame['total_counts'] = np.asarray(totals, dtype=np.int64)
            frame['detected_genes'] = detected
            frame['zero_fraction'] = 1 - frame.detected_genes / len(genes)
            obs.append(frame.reset_index(drop=True))
            inputs.append({'file': path.name, 'sha256': sha, 'cells': len(obj), 'genes': len(genes), 'nonfinite_values': 0})
            offset += len(obj)
            print(f'Ranked context {context}: {offset:,} cells', flush=True)
        finally:
            obj.file.close()
    cells = pd.concat(obs, ignore_index=True)
    if cells.cell_id.duplicated().any():
        raise ValueError('Duplicate cell IDs')
    pd.DataFrame({'gene_name': genes}).to_csv(OUT / 'gene_axis.csv', index=False)
    rng = np.random.default_rng(seed)
    sample = np.sort(rng.choice(len(z), 1500, replace=False))
    sampled = z[sample]
    distance = np.clip(1 - sampled @ sampled.T, 0, 2)
    np.fill_diagonal(distance, 0)
    pd.DataFrame({'cell_id': cells.cell_id.iloc[sample], 'sample_index': sample}).to_csv(OUT / 'silhouette_sample.csv', index=False)
    metrics = []; centers = {}; histories = {}; labels_by_fit = {}
    with threadpool_limits(limits=threads):
        for k, repeat in [(2, 0), (3, 0), (4, 0), (5, 0), (6, 0), (3, 1)]:
            print(f'Fitting full-gene spherical clusters K={k}, repeat={repeat}', flush=True)
            labels, center, history, converged = cluster(z, k, seed + repeat)
            if not converged:
                raise ValueError('Spherical refinement did not converge')
            key = f'k{k}_repeat{repeat}'
            labels_by_fit[key] = labels
            centers[key] = center
            histories[key] = history
            cells['cluster_' + (f'k{k}' if repeat == 0 else key)] = labels
            counts = pd.crosstab(cells.context, labels)
            metric = {'k': k, 'repeat': repeat, 'seed': seed + repeat,
                      'iterations': len(history), 'converged': converged,
                      'sample_silhouette_spearman': separation(sampled, labels[sample], distance),
                      'context_ARI': float(adjusted_rand_score(cells.context, labels)),
                      'context_NMI': float(normalized_mutual_info_score(cells.context, labels)),
                      'context_purity': float(counts.max(axis=0).sum() / len(z)),
                      'min_cluster_cells': int(counts.sum(axis=0).min()),
                      'singleton_clusters': int((counts.sum(axis=0) == 1).sum()),
                      'sample_clusters_represented': len(np.unique(labels[sample])),
                      'sum_centroid_rho': history[-1]['sum_centroid_rho']}
            metrics.append(metric)
            print(metric, flush=True)
        primary = labels_by_fit['k3_repeat0']
        similarity = z @ centers['k3_repeat0'].T
        ordered = np.sort(similarity, axis=1)
        cells['own_cluster_centroid_rho'] = ordered[:, -1]
        cells['cluster_centroid_rho_margin'] = ordered[:, -1] - ordered[:, -2]
        # Exact all-peer group means, with no dense cell-pair matrix.
        pair_means = []
        for target in CONTEXTS:
            member = cells.context.to_numpy() == target
            sums = z[member].sum(axis=0, dtype=np.float64)
            dot = z @ sums
            own = np.einsum('ij,ij->i', z, z, dtype=np.float64)
            mean = (dot - member * own) / (int(member.sum()) - member.astype(int))
            cells[f'mean_peer_rho_{target}'] = mean
            for context in CONTEXTS:
                pair_means.append({'context': context, 'peer_context': target,
                                   'mean_full_gene_spearman': float(mean[cells.context == context].mean())})
        print('Computing label-blind PCA for display only', flush=True)
        pca = PCA(n_components=8, svd_solver='randomized', random_state=seed,
                  iterated_power=3, n_oversamples=10)
        pca.fit(z)
        scores = pca.transform(z)
        for pc in range(8):
            cells[f'PC{pc+1}'] = scores[:, pc]
        centers['pca_components'] = pca.components_
        centers['pca_mean'] = pca.mean_
    metrics = pd.DataFrame(metrics)
    composition = pd.crosstab(cells.context, primary).rename_axis(index='context', columns='cluster').stack().reset_index(name='cells')
    composition['within_context_fraction'] = composition.cells / composition.groupby('context').cells.transform('sum')
    composition['within_cluster_fraction'] = composition.cells / composition.groupby('cluster').cells.transform('sum')
    all_compositions = []
    for key, labels in labels_by_fit.items():
        frame = pd.crosstab(cells.context, labels).rename_axis(index='context', columns='cluster').stack().reset_index(name='cells')
        frame.insert(0, 'fit', key)
        all_compositions.append(frame)
    dispersion = []
    for context in CONTEXTS:
        sub = cells[cells.context == context]
        for name in ['total_counts', 'detected_genes', 'own_cluster_centroid_rho', 'cluster_centroid_rho_margin']:
            values = sub[name].to_numpy()
            dispersion.append({'context': context, 'feature': name, 'mean': float(values.mean()),
                               'sd': float(values.std(ddof=1)), 'p05': float(np.quantile(values, .05)),
                               'median': float(np.median(values)), 'p95': float(np.quantile(values, .95))})
    depth = []; guide = []
    for context in CONTEXTS:
        sub = cells[cells.context == context].copy()
        sub['depth_quintile'] = pd.qcut(sub.total_counts.rank(method='first'), 5, labels=False) + 1
        for quintile in range(1, 6):
            group = sub[sub.depth_quintile == quintile]
            for k in range(3):
                count = int((group.cluster_k3 == k).sum())
                depth.append({'context': context, 'depth_quintile': quintile, 'cluster': k,
                              'cells': count, 'within_depth_fraction': count / len(group)})
        for name, group in sub.groupby('ntc_id', observed=True):
            for k in range(3):
                count = int((group.cluster_k3 == k).sum())
                guide.append({'context': context, 'ntc_id': name, 'cluster': k,
                              'cells': count, 'within_guide_fraction': count / len(group)})
    confounding = []
    for context in CONTEXTS:
        sub = cells[cells.context == context]
        for pc in range(1, 9):
            for feature in ['total_counts', 'detected_genes']:
                confounding.append({'context': context, 'PC': pc, 'feature': feature,
                                    'spearman_rho': float(spearmanr(sub[f'PC{pc}'], sub[feature]).statistic)})
    np.savez_compressed(OUT / 'cluster_models.npz', **centers)
    outputs = {'per_cell_full_rank_clusters.csv': cells, 'cluster_metrics.csv': metrics,
               'cluster_composition_k3.csv': composition,
               'cluster_composition_all_fits.csv': pd.concat(all_compositions, ignore_index=True),
               'within_context_dispersion.csv': pd.DataFrame(dispersion), 'depth_quintile_composition.csv': pd.DataFrame(depth),
               'guide_composition.csv': pd.DataFrame(guide), 'pc_depth_correlations.csv': pd.DataFrame(confounding),
               'all_peer_context_mean_correlations.csv': pd.DataFrame(pair_means),
               'pca_variance.csv': pd.DataFrame({'PC': np.arange(1, 9), 'explained_variance_ratio': pca.explained_variance_ratio_})}
    for name, frame in outputs.items():
        frame.to_csv(OUT / name, index=False)
    report = {'created_utc': datetime.now(timezone.utc).isoformat(), 'seed': seed, 'threads': threads,
              'cells': len(cells), 'genes': len(genes), 'inputs': inputs,
              'method': 'Full-gene unit-centered average ranks, MiniBatchKMeans initialization, exact spherical refinement',
              'metric': 'Spearman rho = dot(unit-centered ranks); silhouette distance = 1 - rho',
              'missing_policy': 'All provided values verified finite; reject any missing value rather than impute',
              'transform': 'No pseudocount, log, CPM, gene filtering, or gene standardization; zeros are tied valid values',
              'labels_used_for_fit': False, 'pca_used_for_fit': False,
              'silhouette_sample_cells': len(sample), 'silhouette_context_labels': separation(sampled, cells.context.iloc[sample], distance),
              'k3_repeat_ARI': float(adjusted_rand_score(primary, labels_by_fit['k3_repeat1'])),
              'best_sample_silhouette_k': int(metrics[metrics.repeat == 0].sort_values(['sample_silhouette_spearman', 'k'], ascending=[False, True]).iloc[0].k),
              'guide_within_context_NMI': {c: float(normalized_mutual_info_score(cells.loc[cells.context == c, 'ntc_id'],
                                                                cells.loc[cells.context == c, 'cluster_k3'])) for c in CONTEXTS},
              'versions': {name: version(name) for name in ['numpy', 'scipy', 'pandas', 'anndata', 'scikit-learn', 'matplotlib']},
              'python': platform.python_version(), 'elapsed_seconds': time.perf_counter() - started,
              'scripts': {name: digest(ROOT / name) for name in ['cluster_full_rank_matrix.py', 'test_full_rank_clustering.py', 'validate_full_rank_clustering.py']},
              'refinement_history': histories}
    make_plot(cells, composition, metrics, pca.explained_variance_ratio_)
    report['outputs'] = {p.name: {'bytes': p.stat().st_size, 'sha256': digest(p)} for p in OUT.iterdir()
                         if p.is_file() and p.name not in {'analysis_report.json', 'validation.json'}}
    (OUT / 'analysis_report.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    write_results(metrics, composition, report, pair_means, confounding, dispersion)
    print('Completed full rank clustering', flush=True)


def write_results(metrics, composition, report, pair_means, confounding, dispersion):
    text = '# Complete-matrix rank clustering results\n\n'
    text += 'All 55,200 supplied non-targeting cells and all 18,533 measured genes were clustered together. A/B/C labels were withheld from fitting. Cluster numbers are arbitrary.\n\n'
    text += '## Three-cluster composition\n\n' + table(composition) + '\n\n'
    if (metrics.loc[metrics.k == 3, 'context_ARI'] == 1).all():
        text += 'Both independent K=3 fits assign every A cell to one pure cluster, every B cell to another, and every C cell to a third: **no context mixing among any of the 55,200 cells**. ARI, NMI and purity are all 1.0.\n\n'
    text += '## Number-of-clusters sensitivity\n\n' + table(metrics[['k', 'repeat', 'sample_silhouette_spearman', 'context_ARI', 'context_purity', 'min_cluster_cells', 'singleton_clusters', 'sample_clusters_represented']]) + '\n\n'
    text += f"Three-cluster repeat ARI: **{report['k3_repeat_ARI']:.6f}**. The largest sampled silhouette among K=2..6 first occurred at K={report['best_sample_silhouette_k']} (ties resolved toward smaller K). The silhouette of the supplied context labels on the same sample was {report['silhouette_context_labels']:.6f}. This finite candidate range does not establish the true number of biological cell states.\n\n"
    text += 'Higher-K singleton counts and the number of clusters represented in the silhouette sample are shown above. Omitted singleton clusters can make different full partitions receive identical sampled silhouettes. The complete nonzero composition is below; higher-K splits and local optima need further study before biological interpretation.\n\n'
    all_fits = pd.read_csv(OUT / 'cluster_composition_all_fits.csv')
    text += table(all_fits[all_fits.cells > 0]) + '\n\n'
    text += '## Within-context dispersion\n\n' + table(pd.DataFrame(dispersion)) + '\n\n'
    text += 'Count-depth distributions and cluster composition within depth quintiles and guides are recorded. When all cells have context-pure assignments, that purity necessarily persists in each stratum. The projection also shows depth gradients within contexts; neither clustering nor similar count depths removes dropout or establishes exclusively biological separation.\n\n'
    text += '## Exact mean correlations to all context peers\n\n' + table(pd.DataFrame(pair_means)) + '\n\n'
    text += 'Every peer is included; self-pairs are excluded. This calculation uses full-gene dot products of unit rank vectors and group sums, without sampling or a dense 55,200-square matrix.\n\n'
    text += '## Depth and guide diagnostics\n\n' + table(pd.DataFrame(confounding).query('PC <= 2')) + '\n\n'
    text += 'Within-context guide/cluster NMI: ' + ', '.join(f"{c}: {v:.6f}" for c,v in report['guide_within_context_NMI'].items()) + '. See depth-quintile and guide composition CSVs for all cells. PC axes describe only a projection, and count-depth correlations do not establish causation.\n\n'
    text += 'These data contain only supplied non-targeting controls. Context separation can include biological and batch differences; it cannot identify the anonymous cell types or validate CRISPRi responses.\n\n'
    text += '![Full rank clustering](analysis/full_rank_clustering/full_rank_clusters.png)\n\nReproduction, methods and validation: [FULL_RANK_CLUSTERING.md](FULL_RANK_CLUSTERING.md).\n'
    (ROOT / 'FULL_RANK_CLUSTERING_RESULTS.md').write_text(text, encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, default=20261004)
    parser.add_argument('--threads', type=int, default=12)
    args = parser.parse_args()
    run(args.seed, args.threads)
