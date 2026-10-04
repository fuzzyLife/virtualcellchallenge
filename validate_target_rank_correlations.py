"""Independently check target-only rank results against original counts and SciPy."""
from pathlib import Path
import hashlib
import json
import warnings
import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, ConstantInputWarning

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'analysis/target_rank_correlations'
REFERENCES = ['K562', 'HCT116', 'HEK293T', 'H1', 'CD4']


def oracle(left, right, minimum):
    valid = np.isfinite(left) & np.isfinite(right)
    count = int(valid.sum())
    if count < minimum:
        return np.nan, count
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', ConstantInputWarning)
        return float(spearmanr(left[valid], right[valid]).statistic), count


def equal(left, right):
    return bool(np.isclose(left, right, atol=1e-12, rtol=0, equal_nan=True))


def main():
    report = json.loads((OUT / 'analysis_report.json').read_text())
    cells = pd.read_csv(OUT / 'per_cell_target_rank_correlations.csv')
    summary = pd.read_csv(OUT / 'context_reference_summary.csv')
    targets = pd.read_csv(ROOT / 'pert_counts.csv').target_gene.to_numpy()
    genes = pd.read_csv(ROOT / 'gene_names.csv').gene_name.to_numpy()
    positions = {gene: i for i, gene in enumerate(genes)}
    axis = np.array([positions[gene] for gene in targets])
    reference = pd.read_csv(OUT / 'target_reference_raw_profiles.csv')
    assert np.array_equal(reference.gene_name.to_numpy(), targets)
    assert len(targets) == 300
    for name in REFERENCES:
        prepared = pd.read_csv(ROOT / f'reference/control_profiles/{name}.csv').set_index('gene_name').control_count_sum.reindex(targets).to_numpy()
        np.testing.assert_allclose(reference[name].to_numpy(), prepared, rtol=1e-14, atol=0, equal_nan=True)
    checks = []
    for context in ['A', 'B', 'C']:
        obj = ad.read_h5ad(ROOT / f'context_{context}.h5ad', backed='r')
        sub = cells[cells.context == context]
        assert len(sub) == obj.n_obs == 18400
        assert np.array_equal(sub.cell_id.to_numpy(), obj.obs_names.to_numpy())
        assert np.array_equal(sub.ntc_id.to_numpy(), obj.obs.ntc_id.astype(str).to_numpy())
        for name in REFERENCES:
            values = sub['spearman_' + name]
            finite = values.notna()
            assert (sub.loc[finite, 'status_' + name] == 'ok').all()
            assert (sub.loc[~finite, 'status_' + name] != 'ok').all()
            assert values.dropna().between(-1, 1).all()
            expected_count = int(np.isfinite(reference[name]).sum())
            assert (sub['valid_gene_pairs_' + name] == expected_count).all()
            entry = summary[(summary.context == context) & (summary.reference == name)].iloc[0]
            assert equal(values.median(), entry.median_spearman)
            assert int(finite.sum()) == entry.defined_correlations
        for index in [0, 9199, 18399]:
            raw = obj.X[index:index + 1].toarray().ravel().astype(np.float64)
            target_counts = raw[axis]
            entry = sub.iloc[index]
            assert entry.total_counts == raw.sum()
            assert entry.nonzero_target_genes == np.count_nonzero(target_counts)
            for name in REFERENCES:
                expected, count = oracle(target_counts, reference[name].to_numpy(), report['min_valid_pairs'])
                assert equal(expected, entry['spearman_' + name]), (context, index, name)
                assert count == entry['valid_gene_pairs_' + name]
                checks.append({'context': context, 'cell_index_0based': index, 'reference': name, 'spearman': expected if np.isfinite(expected) else None, 'valid_gene_pairs': count})
        obj.file.close()
    pooled = pd.read_csv(OUT / 'pooled_raw_target_profiles.csv')
    matrix = pd.read_csv(OUT / 'pooled_profile_rank_correlations.csv', index_col='profile')
    counts = pd.read_csv(OUT / 'pooled_profile_valid_gene_counts.csv', index_col='profile')
    pair_checks = 0
    for left in matrix.index:
        for right in matrix.columns:
            expected, count = oracle(pooled[left].to_numpy(), pooled[right].to_numpy(), report['min_valid_pairs'])
            assert equal(expected, matrix.loc[left, right]), (left, right)
            assert count == counts.loc[left, right]
            pair_checks += 1
    assert len(cells) == 55200 and not cells[['context', 'cell_id']].duplicated().any()
    for name, record in report['outputs'].items():
        assert hashlib.sha256((OUT / name).read_bytes()).hexdigest() == record['sha256'], name
    result = {'passed': True, 'complete_rows': len(cells), 'all_cell_ids_and_guides_match': True, 'all_summary_medians_and_defined_counts_match': True, 'all_pairwise_valid_gene_counts_match': True, 'all_output_csv_checksums_match': True, 'independent_actual_cell_pairs_checked': len(checks), 'pooled_matrix_pairs_checked': pair_checks, 'actual_cell_checks': checks, 'validator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT / 'validation.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    print('Validated all 55,200 rows, 45 original-count Spearman checks and 64 pooled-profile pairs.')


if __name__ == '__main__':
    main()
