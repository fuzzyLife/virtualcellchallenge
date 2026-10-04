"""Independently verify all archive counts, identities, model rows and marker tables."""
from pathlib import Path
import json
import hashlib
import zipfile
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.stats import rankdata
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, roc_auc_score

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'analysis/target_archive_audit'
def sha(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024**2),b''):result.update(block)
    return result.hexdigest()

def main():
    report=json.loads((OUT/'analysis_report.json').read_text())
    audit=json.loads((OUT/'archive_validation.json').read_text())
    assert sha(ROOT/audit['file'])==report['archive_sha256']==audit['sha256']
    for name,value in report['scripts'].items():assert sha(ROOT/name)==value,name
    for name,value in report['outputs'].items():assert sha(OUT/name)==value['sha256'],name
    targets=pd.read_csv(ROOT/'pert_counts.csv').target_gene.to_numpy()
    genes=pd.read_csv(ROOT/'gene_names.csv').gene_name.to_numpy()
    cells=pd.read_csv(OUT/'per_cell_target300_clusters.csv')
    broad=pd.read_csv(OUT/'broader_lineage_expression_summary.csv')
    bg=broad.gene_name.drop_duplicates().to_numpy()
    positions={g:i for i,g in enumerate(genes)}
    with zipfile.ZipFile(ROOT/audit['file']) as archive:
        assert archive.testzip() is None
        assert hashlib.sha256(archive.read(audit['member'])).hexdigest()==audit['member_sha256']
        with archive.open(audit['member']) as handle:frame=pd.read_csv(handle)
    assert list(frame.columns)==['context']+list(targets)
    assert frame.context.tolist()==cells.context.tolist()
    assert np.array_equal(pd.read_csv(OUT/'gene_axis.csv').gene_name,targets)
    raw=frame[targets].to_numpy(float)
    assert raw.shape==(55200,300) and np.isfinite(raw).all() and (raw>=0).all() and (raw==np.floor(raw)).all()
    libraries=np.empty(len(cells));braw=np.empty((len(cells),len(bg)));checked=0
    for record in audit['original_inputs']:
        path=ROOT/record['file'];assert sha(path)==record['sha256']
        obj=ad.read_h5ad(path,backed='r');idx=np.flatnonzero(cells.context.to_numpy()==record['context'])
        try:
            assert obj.shape==(18400,18533) and np.array_equal(obj.var_names,genes)
            assert np.array_equal(obj.obs_names,cells.loc[idx,'cell_id'])
            assert np.array_equal(obj.obs.ntc_id.astype(str),cells.loc[idx,'ntc_id'].astype(str))
            for start in range(0,len(obj),512):
                end=min(start+512,len(obj));x=sparse.csr_matrix(obj.X[start:end],dtype=float);rows=idx[start:end]
                expected=x[:,[positions[g] for g in targets]].toarray()
                assert np.array_equal(expected,raw[rows]),f'{path.name} rows {start}:{end}'
                libraries[rows]=np.asarray(x.sum(axis=1)).ravel()
                braw[rows]=x[:,[positions[g] for g in bg]].toarray();checked+=expected.size
        finally:obj.file.close()
        print(f'Independent complete count comparison passed: {record["context"]}',flush=True)
    assert checked==audit['values_checked']==16560000 and audit['mismatched_values']==0
    assert np.array_equal(libraries,cells.full_library_counts)
    assert np.array_equal(raw.sum(axis=1),cells.target_count_sum)
    assert np.array_equal((raw>0).sum(axis=1),cells.nonzero_target_genes)
    assert cells.cell_id.is_unique
    models=np.load(OUT/'target300_rank_models.npz')
    sample=np.unique(np.concatenate([np.random.default_rng(1701).choice(len(cells),60,replace=False),np.flatnonzero(~cells.matches_context)]))
    ranks=rankdata(raw[sample],method='average',axis=1);ranks-=ranks.mean(axis=1,keepdims=True);ranks/=np.linalg.norm(ranks,axis=1,keepdims=True)
    for key in ['k2_repeat0','k3_repeat0','k4_repeat0','k6_repeat0','k3_repeat1']:
        assert np.array_equal((ranks@models[key].T).argmax(axis=1),cells.loc[sample,'cluster_'+key])
    pcs=(ranks-models['pca_mean'])@models['pca_components'].T
    np.testing.assert_allclose(pcs,cells.loc[sample,[f'PC{i}' for i in range(1,9)]],atol=2e-6,rtol=1e-5)
    sim=ranks@models['context_order_centroids'].T
    assert np.array_equal(np.array(['A','B','C'])[sim.argmax(axis=1)],cells.loc[sample,'assigned_context'])
    np.testing.assert_allclose(sim.max(axis=1),cells.loc[sample,'best_centroid_similarity'],atol=2e-6)
    np.testing.assert_allclose(np.sort(sim,axis=1)[:,-1]-np.sort(sim,axis=1)[:,-2],cells.loc[sample,'centroid_margin'],atol=2e-6)
    matches=cells.context==cells.assigned_context
    assert np.array_equal(matches,cells.matches_context)
    assert matches.sum()==report['correctly_grouped_cells']
    assert (~matches).sum()==report['disagreeing_cells']
    composition=pd.read_csv(OUT/'cluster_composition.csv')
    for row in composition.itertuples():assert row.cells==((cells.context==row.context)&(cells.assigned_context==row.assigned_context)).sum()
    metrics=pd.read_csv(OUT/'cluster_metrics.csv')
    for row in metrics.itertuples():
        labels=cells[f'cluster_k{row.k}_repeat{row.repeat}']
        np.testing.assert_allclose(row.context_ARI,adjusted_rand_score(cells.context,labels),atol=1e-14)
        np.testing.assert_allclose(row.context_NMI,normalized_mutual_info_score(cells.context,labels),atol=1e-14)
    np.testing.assert_allclose(report['repeat_ARI'],adjusted_rand_score(cells.cluster_k3_repeat0,cells.cluster_k3_repeat1),atol=1e-14)
    marker_checks=0
    for filename,column in [('target300_cluster_markers.csv','assigned_context'),('target300_context_markers.csv','context')]:
        markers=pd.read_csv(OUT/filename);groups=cells[column].to_numpy();cpm=raw/libraries[:,None]*1e6
        assert len(markers)==900 and not markers[['group','gene_name']].duplicated().any()
        assert set(markers.gene_name)==set(targets)
        for row in markers.itertuples():
            j=list(targets).index(row.gene_name);mask=groups==row.group;other=[g for g in ['A','B','C'] if g!=row.group]
            mean=cpm[mask,j].mean();den=max(cpm[groups==g,j].mean() for g in other)
            expected_fold=mean/den if den>0 else (np.inf if mean>0 else np.nan)
            aucs=[roc_auc_score(np.r_[np.ones(mask.sum()),np.zeros((groups==g).sum())],np.r_[cpm[mask,j],cpm[groups==g,j]]) for g in other]
            np.testing.assert_allclose([row.mean_raw_count,row.median_raw_count,row.detection_fraction,row.mean_cpm,row.median_cpm,row.max_other_mean_cpm,row.fold_over_highest_other_group,row.minimum_pairwise_auc,row.auc_vs_rest],
               [raw[mask,j].mean(),np.median(raw[mask,j]),(raw[mask,j]>0).mean(),mean,np.median(cpm[mask,j]),den,expected_fold,min(aucs),roc_auc_score(mask,cpm[:,j])],rtol=1e-10,atol=1e-10,equal_nan=True)
            assert row.panel_enriched==(min(aucs)>=.8 and expected_fold>=2 and (raw[mask,j]>0).mean()>=.25)
            assert row.detection_restricted==((raw[mask,j]>0).mean()>=.5 and max((raw[groups==g,j]>0).mean() for g in other)<=.05)
            marker_checks+=1
    for row in broad.itertuples():
        mask=cells.context.to_numpy()==row.context;j=list(bg).index(row.gene_name);x=braw[mask,j];cpm=x/libraries[mask]*1e6
        np.testing.assert_allclose([row.mean_raw_count,row.median_raw_count,row.detection_fraction,row.mean_cpm,row.median_cpm],[x.mean(),np.median(x),(x>0).mean(),cpm.mean(),np.median(cpm)],rtol=1e-10,atol=1e-10)
    discrep=pd.read_csv(OUT/'cells_disagreeing_with_context.csv')
    assert set(discrep.cell_id)==set(cells.loc[~matches,'cell_id'])
    original=pd.read_csv(ROOT/'analysis/full_rank_clustering/per_cell_full_rank_clusters.csv')
    for row in discrep.itertuples():
        source=original[original.cell_id==row.cell_id].iloc[0]
        assert row.original_full_gene_context==source.context and row.original_full_gene_cluster==source.cluster_k3
    result={'status':'passed','all_original_count_values_checked':checked,'all_cell_identities_and_depths_checked':len(cells),
            'all_gene_group_marker_rows_independently_checked':marker_checks,'broader_marker_rows_checked':len(broad),
            'original_count_rank_and_model_cells_checked':len(sample),'output_and_script_checksums_verified':True}
    (OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
