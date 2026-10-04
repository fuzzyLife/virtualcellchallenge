"""Independent checks of cell dispersion using original counts and SciPy ranks."""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.stats import rankdata, spearmanr

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'analysis/target_cell_dispersion'
REFS=['K562','HCT116','HEK293T','H1','CD4']


def main():
    report=json.loads((OUT/'analysis_report.json').read_text())
    cells=pd.read_csv(OUT/'per_cell_dispersion.csv')
    original=pd.read_csv(ROOT/'analysis/target_rank_correlations/per_cell_target_rank_correlations.csv')
    assert len(cells)==55200 and not cells[['context','cell_id']].duplicated().any()
    pd.testing.assert_frame_equal(cells[original.columns],original)
    genes=pd.read_csv(ROOT/'gene_names.csv').gene_name.to_numpy()
    targets=pd.read_csv(ROOT/'pert_counts.csv').target_gene.to_numpy()
    pos={g:i for i,g in enumerate(genes)}; axis=[pos[g] for g in targets]
    reference=pd.read_csv(ROOT/'analysis/target_rank_correlations/target_reference_raw_profiles.csv')
    common=np.isfinite(reference[REFS].to_numpy()).all(axis=1)
    gene_axis=pd.read_csv(OUT/'gene_axis.csv')
    np.testing.assert_array_equal(gene_axis.gene_name,targets)
    np.testing.assert_array_equal(gene_axis.in_common_reference_axis,common)
    assert int(common.sum())==231
    raw_arrays=[]
    gene_stats=pd.read_csv(OUT/'target_gene_count_dispersion.csv')
    for context in ['A','B','C']:
        obj=ad.read_h5ad(ROOT/f'context_{context}.h5ad',backed='r')
        sub=cells[cells.context==context]
        np.testing.assert_array_equal(sub.cell_id,obj.obs_names)
        np.testing.assert_array_equal(sub.ntc_id,obj.obs.ntc_id.astype(str))
        raw=np.empty((obj.n_obs,300))
        for start in range(0,obj.n_obs,1024):
            end=min(start+1024,obj.n_obs)
            raw[start:end]=sparse.csr_matrix(obj.X[start:end],dtype=float)[:,axis].toarray()
        raw_arrays.append(raw); obj.file.close()
        stats=gene_stats[gene_stats.context==context]
        np.testing.assert_array_equal(stats.gene_name,targets)
        np.testing.assert_allclose(stats.mean_raw_count,raw.mean(axis=0),rtol=1e-12)
        np.testing.assert_allclose(stats.sample_variance_raw_count,raw.var(axis=0,ddof=1),rtol=1e-12)
        np.testing.assert_allclose(stats.zero_fraction,(raw==0).mean(axis=0),atol=1e-14)
        mean=raw.mean(axis=0); var=raw.var(axis=0,ddof=1)
        np.testing.assert_allclose(stats.fano_raw_counts,np.divide(var,mean,out=np.full(300,np.nan),where=mean>0),rtol=1e-12,equal_nan=True)
        print(context,'original counts and variance validated',flush=True)
    raw=np.vstack(raw_arrays)
    ranks=rankdata(raw,axis=1,method='average')
    ranks-=ranks.mean(axis=1,keepdims=True)
    norms=np.linalg.norm(ranks,axis=1)
    assert (norms>0).all() and report['constant_cell_vectors']==0
    ranks/=norms[:,None]
    peer_checks=0; profile_checks=0; loo_checks=0
    for ci,context in enumerate(['A','B','C']):
        indexes=ci*18400+np.array([0,9199,18399])
        for i in indexes:
            row=cells.iloc[i]
            for peer in ['A','B','C']:
                ids=np.flatnonzero(cells.context.to_numpy()==peer)
                ids=ids[ids!=i]
                # Explicitly construct every correlation to this cell, rather than using moments.
                all_rho=ranks[ids]@ranks[i]
                np.testing.assert_allclose(row['mean_cell_correlation_to_'+peer],all_rho.mean(),atol=1e-12,rtol=0)
                np.testing.assert_allclose(row['sd_cell_correlation_to_'+peer],all_rho.std(ddof=0),atol=1e-12,rtol=0)
                assert row['valid_peer_cells_in_'+peer]==len(ids)
                peer_checks+=1
            for ref in REFS:
                expected=spearmanr(raw[i,common],reference.loc[common,ref]).statistic
                np.testing.assert_allclose(row['common231_spearman_'+ref],expected,atol=1e-12,rtol=0)
                profile_checks+=1
            own=raw_arrays[ci].sum(axis=0)-raw[i]
            expected=spearmanr(raw[i],own).statistic
            np.testing.assert_allclose(row.spearman_own_context_leave_one_out_300,expected,atol=1e-12,rtol=0)
            loo_checks+=1
    for peer in ['A','B','C']:
        np.testing.assert_array_equal(cells['valid_peer_cells_in_'+peer],np.where(cells.context==peer,18399,18400))
        assert cells['mean_cell_correlation_to_'+peer].between(-1,1).all()
        assert cells['sd_cell_correlation_to_'+peer].between(0,1).all()
    np.testing.assert_allclose(cells.delta_CD4_minus_H1_300,cells.spearman_CD4-cells.spearman_H1,atol=1e-14)
    expected_winner=np.array(REFS)[cells[['common231_spearman_'+ref for ref in REFS]].to_numpy().argmax(axis=1)]
    np.testing.assert_array_equal(cells.highest_common231_reference,expected_winner)
    pairs=pd.read_csv(OUT/'sampled_individual_cell_pair_correlations.csv')
    assert len(pairs)==6*report['sample_pairs_per_group']
    mapping={cell:i for i,cell in enumerate(cells.cell_id)}
    assert len(mapping)==len(cells)
    left=np.array([mapping[x] for x in pairs.cell_id_left]); right=np.array([mapping[x] for x in pairs.cell_id_right])
    assert (left!=right).all() and (pairs.valid_gene_pairs==300).all()
    np.testing.assert_array_equal(pairs.cell_context_left,cells.context.to_numpy()[left])
    np.testing.assert_array_equal(pairs.cell_context_right,cells.context.to_numpy()[right])
    expected=np.einsum('ij,ij->i',ranks[left],ranks[right])
    np.testing.assert_allclose(pairs.spearman,expected,atol=1e-12,rtol=0)
    for i in np.linspace(0,len(pairs)-1,36,dtype=int):
        np.testing.assert_allclose(pairs.iloc[i].spearman,spearmanr(raw[left[i]],raw[right[i]]).statistic,atol=1e-12,rtol=0)
    summary=pd.read_csv(OUT/'correlation_dispersion_summary.csv')
    for entry in summary.itertuples():
        sub=cells[cells.context==entry.context]
        if entry.gene_axis=='pairwise_available_targets':
            col='spearman_'+entry.comparison
        elif entry.gene_axis=='common_231_targets':
            col='common231_spearman_'+entry.comparison
        elif entry.comparison.startswith('mean_to_all_cells_'):
            col='mean_cell_correlation_to_'+entry.comparison[-1]
        else:
            col=entry.comparison
        x=sub[col].dropna().to_numpy()
        np.testing.assert_allclose([entry.mean,entry.variance,entry.sd,entry.p05,entry.median,entry.p95],
            [x.mean(),x.var(ddof=1),x.std(ddof=1),*np.quantile(x,[.05,.5,.95])],atol=1e-12,rtol=0)
    pair_summary=pd.read_csv(OUT/'cell_pair_dispersion_summary.csv')
    for entry in pair_summary.itertuples():
        a,b=entry.pair_group.split('-')
        sub=cells[cells.context==a]
        count=sub['valid_peer_cells_in_'+b].to_numpy()
        mean=sub['mean_cell_correlation_to_'+b].to_numpy()
        sd=sub['sd_cell_correlation_to_'+b].to_numpy()
        m=np.average(mean,weights=count)
        s=np.sqrt(np.average(sd**2+mean**2,weights=count)-m*m)
        np.testing.assert_allclose([entry.all_pairs_mean,entry.all_pairs_population_sd],[m,s],atol=1e-12,rtol=0)
        assert entry.all_unique_valid_cell_pairs==count.sum()//(2 if a==b else 1)
    for name,record in report['outputs'].items():
        assert hashlib.sha256((OUT/name).read_bytes()).hexdigest()==record['sha256'],name
    result={'passed':True,'all_cell_rows_checked':len(cells),'sampled_pair_rows_checked':len(pairs),
        'explicit_all_peer_mean_sd_checks':peer_checks,'independent_shared_axis_scipy_checks':profile_checks,
        'independent_leave_one_out_scipy_checks':loo_checks,'individual_pair_scipy_checks':36,
        'target_gene_variance_rows_checked':len(gene_stats),'output_checksums_verified':True,
        'validator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT/'validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    main()
