"""Small response kernels with pairwise-valid rank context features.

A row is a measured target/context contrast. Output masks distinguish unavailable
responses from measured zero. All panel genes are excluded from fitted outputs.
"""
import numpy as np
from scipy import linalg
from scipy.stats import rankdata
from analyze_target_rank_correlations import pairwise_spearman

VARIANTS = ('identity_rank', 'pca_rank', 'tahoe_rank', 'tahoe_pool')
RIDGE = 10.0


def normalize_rows(x):
    x = np.asarray(x, dtype=float)
    norm = np.linalg.norm(x, axis=1, keepdims=True)
    return np.divide(x, norm, out=np.zeros_like(x), where=norm > 1e-12)


def rank_context_features(profiles, training_contexts):
    refs = np.asarray(training_contexts, dtype=int)
    result = np.zeros((len(profiles), len(refs)))
    counts = np.zeros_like(result, dtype=int)
    for i, profile in enumerate(profiles):
        for j, ref in enumerate(refs):
            rho, n, status = pairwise_spearman(profile, profiles[ref], min_pairs=100)
            if status != 'ok':
                raise ValueError('Insufficient nonconstant control profile for rank kernel')
            result[i, j], counts[i, j] = rho, n
    return result, counts


def pooled_context_features(profiles, gene_vectors):
    # Engineered pooling, not inference through a released Tahoe cell encoder.
    vectors, available = gene_vectors
    result = []
    for profile in profiles:
        valid = available & np.isfinite(profile)
        if valid.sum() < 100:
            raise ValueError('Insufficient pooling genes')
        ranks = rankdata(profile[valid], method='average')
        ranks -= ranks.mean()
        ranks /= np.linalg.norm(ranks)
        result.append(ranks @ vectors[valid])
    return normalize_rows(np.asarray(result))


def kernels(row_targets, row_contexts, query_targets, query_contexts,
            profiles, embeddings, variant, gene_vectors=None):
    row_targets, query_targets = np.asarray(row_targets), np.asarray(query_targets)
    refs = np.unique(row_contexts)
    if variant == 'tahoe_pool':
        features = pooled_context_features(profiles, gene_vectors)
    else:
        features, _ = rank_context_features(profiles, refs)
    # Center only on distinct training contexts. No treated profiles enter here.
    features = features - features[refs].mean(axis=0)
    train_c, query_c = features[row_contexts], features[query_contexts]
    def rbf(a, b):
        distance = np.maximum((a*a).sum(axis=1)[:, None] + (b*b).sum(axis=1)[None, :] - 2*a@b.T, 0)
        return np.exp(-2*distance)
    context_train = .25 + .75*rbf(train_c, train_c)
    context_query = .25 + .75*rbf(query_c, train_c)
    target_train = (row_targets[:, None] == row_targets[None, :]).astype(float) + .1
    target_query = (query_targets[:, None] == row_targets[None, :]).astype(float) + .1
    if variant != 'identity_rank':
        family = 'pca' if variant == 'pca_rank' else 'tahoe'
        vectors, available = embeddings[family]
        unique = np.unique(row_targets)
        usable = unique[available[unique]]
        if len(usable) < 5:
            raise ValueError('Too few training embeddings')
        centre = vectors[usable].mean(axis=0)
        scale = vectors[usable].std(axis=0)
        scale = np.where(scale > 1e-8, scale, 1.)
        transformed = normalize_rows((vectors-centre)/scale)
        transformed[~available] = 0  # Explicitly no embedding kernel for missing genes.
        target_train += transformed[row_targets] @ transformed[row_targets].T
        target_query += transformed[query_targets] @ transformed[row_targets].T
    return target_train*context_train, target_query*context_query


def predict_ridge(kernel, cross_kernel, values, valid, row_families, ridge=RIDGE):
    """Solve independently for each output observation mask, never impute labels.

    Each training lineage has equal total weight within each output's valid rows.
    Fixed ridge=10 is not selected on held-out outcomes.
    """
    valid = np.asarray(valid, bool)
    patterns, inverse = np.unique(np.packbits(valid.T, axis=1), axis=0, return_inverse=True)
    prediction = np.zeros((len(cross_kernel), values.shape[1]), np.float32)
    supported = np.zeros(values.shape[1], bool)
    for group in range(len(patterns)):
        columns = np.flatnonzero(inverse == group)
        rows = np.flatnonzero(valid[:, columns[0]])
        if len(rows) < 5:
            continue
        families = np.asarray(row_families)[rows]
        labels, counts = np.unique(families, return_counts=True)
        weight = np.array([len(rows)/len(labels)/counts[np.flatnonzero(labels==f)[0]] for f in families])
        root = np.sqrt(weight)
        system = kernel[np.ix_(rows, rows)]*root[:, None]*root[None, :]
        system.flat[::len(rows)+1] += ridge
        rhs = values[np.ix_(rows, columns)]*root[:, None]
        factor = linalg.cho_factor(system, check_finite=False)
        coefficient = linalg.cho_solve(factor, rhs, check_finite=False)
        prediction[:, columns] = (cross_kernel[:, rows]*root[None, :]) @ coefficient
        supported[columns] = True
    return np.clip(prediction, -3, 3), supported
