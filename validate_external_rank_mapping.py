"""Independent checks of source axes, raw ranks, frozen assignments and summaries."""
import json
import numpy as np
import pandas as pd
import anndata as ad
import h5py
import pyarrow.parquet as pq
from scipy import sparse
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import adjusted_rand_score
from analyze_target_rank_correlations import ROOT, digest, REFERENCES
from extract_external_mapping_controls import CACHE, PROFILE_DIR, CachedRangeReader
from prepare_reference_controls import column

OUT=ROOT/'analysis/external_rank_mapping'
CONTEXTS=['A','B','C']


def independent_projection(raw, model):
    z=rankdata(np.asarray(raw,dtype=np.float64),axis=1,method='average')
    z-=z.mean(axis=1,keepdims=True);z/=np.linalg.norm(z,axis=1,keepdims=True)
    scores=z@model['centroids'].T
    pcs=(z-model['pca_mean'])@model['pca_components'].T
    return z,scores,pcs


def verify_projected(raw, saved, model):
    z, scores, pcs=independent_projection(raw,model)
    np.testing.assert_array_equal(np.asarray(CONTEXTS)[scores.argmax(axis=1)],saved.nearest_context)
    np.testing.assert_allclose(scores,saved[['rank_centroid_similarity_'+c for c in CONTEXTS]],atol=2e-6)
    np.testing.assert_allclose(pcs,saved[[f'PC{i}' for i in range(1,9)]],atol=3e-6)
    if 'mean_cell_spearman_A' in saved:
        np.testing.assert_allclose(z@model['context_mean_ranks'].T,saved[['mean_cell_spearman_'+c for c in CONTEXTS]],atol=2e-6)
        ordered=np.sort(scores,axis=1)
        np.testing.assert_allclose(ordered[:,-1]-ordered[:,-2],saved.centroid_similarity_margin,atol=2e-6)
        distances=[]
        for i,pc in enumerate(pcs):
            c=int(scores[i].argmax());delta=pc-model['pc_context_means'][c]
            distances.append(delta@model['pc_context_cov_inverse'][c]@delta)
        np.testing.assert_allclose(distances,saved.pc8_mahalanobis_squared,rtol=2e-4,atol=2e-4)
        all_distances=[]
        for c in range(3):
            delta=pcs-model['pc_context_means'][c]
            distance=np.einsum('ij,ij->i',delta@model['pc_context_cov_inverse'][c],delta)
            np.testing.assert_allclose(distance,saved['pc8_mahalanobis_squared_'+CONTEXTS[c]],rtol=2e-4,atol=2e-4)
            all_distances.append(distance)
        np.testing.assert_array_equal((np.array(all_distances).T>model['pc_context_mahal_p95']).all(axis=1),saved.outside_all_context_pc8_95pct)
    return z


def run():
    report=json.loads((OUT/'analysis_report.json').read_text())
    extraction=json.loads((OUT/'external_control_extraction.json').read_text())
    for name,record in report['outputs'].items():assert digest(OUT/name)==record['sha256'],name
    for name,sha in report['scripts'].items():assert digest(ROOT/name)==sha,name
    for name,sha in extraction['outputs'].items():assert digest(CACHE/name)==sha,name
    assert digest(ROOT/'analysis/full_rank_clustering/cluster_models.npz')==report['original_full_rank_model_sha256']
    axis=pd.read_csv(OUT/'gene_axis.csv').gene_name.to_numpy()
    genes=pd.read_csv(ROOT/'gene_names.csv').gene_name.to_numpy()
    profiles={n:pd.read_csv(PROFILE_DIR/f'{n}.csv').set_index('gene_name').control_count_sum for n in REFERENCES}
    expected=np.array([g for g in genes if all(g in profiles[n].index for n in REFERENCES)])
    np.testing.assert_array_equal(axis,expected);assert len(axis)==7447
    cells=pd.read_csv(OUT/'competition_cells_shared_space.csv')
    external=pd.read_csv(OUT/'external_control_units_mapped.csv')
    pooled=pd.read_csv(OUT/'pooled_profiles_mapped.csv')
    summary=pd.read_csv(OUT/'external_mapping_summary.csv')
    model=np.load(OUT/'frozen_shared_rank_model.npz',allow_pickle=False)
    raw_cache=sparse.load_npz(CACHE/'counts.npz')
    meta=pd.read_csv(CACHE/'units.csv')
    assert len(cells)==55200 and len(external)==4393 and not external.unit_id.duplicated().any()
    np.testing.assert_array_equal(meta.unit_id,external.unit_id)
    assert abs(adjusted_rand_score(cells.context,cells.blind_cluster)-report['competition_context_ARI'])<1e-12
    count_table=pd.crosstab(cells.context,cells.blind_cluster)
    np.testing.assert_array_equal(count_table,pd.read_csv(OUT/'competition_cluster_composition.csv').pivot(index='context',columns='blind_cluster',values='cells'))
    for row in summary.itertuples():
        sub=external[external.reference==row.reference]
        if row.stratum!='all':sub=sub[sub.culture_condition==row.stratum]
        assert len(sub)==row.units and int(sub.n_cells.sum())==row.represented_cells
        for c in CONTEXTS:assert abs(getattr(row,'fraction_'+c)-(sub.nearest_context==c).mean())<1e-12
        assert abs(row.outside_pc8_fraction-sub.outside_assigned_context_pc8_95pct.mean())<1e-12
        assert abs(row.outside_all_pc8_fraction-sub.outside_all_context_pc8_95pct.mean())<1e-12
        assert abs(row.median_similarity_margin-sub.centroid_similarity_margin.median())<2e-7
    rng=np.random.default_rng(701)
    source_pairs=[]
    for context in CONTEXTS:
        record=next(x for x in report['competition_inputs'] if x['file']==f'context_{context}.h5ad')
        assert digest(ROOT/record['file'])==record['sha256']
        obj=ad.read_h5ad(ROOT/record['file'],backed='r')
        try:
            sub=cells[cells.context==context].reset_index(drop=True)
            np.testing.assert_array_equal(sub.cell_id,obj.obs_names)
            np.testing.assert_array_equal(sub.ntc_id,obj.obs.ntc_id.astype(str))
            positions=obj.var_names.get_indexer(axis);assert (positions>=0).all()
            indices=np.sort(rng.choice(len(obj),10,replace=False))
            raw=obj.X[indices][:,positions].toarray()
            np.testing.assert_allclose(sub.common_total_counts.iloc[indices],raw.sum(axis=1))
            verify_projected(raw,sub.iloc[indices],model)
            source_pairs.extend(raw[:3])
        finally:obj.file.close()
    for name in REFERENCES:
        sub=external[external.reference==name]
        chosen=pd.read_csv(PROFILE_DIR/f'{name}_selected_controls.csv')
        np.testing.assert_array_equal(sub.source_row,chosen.source_row)
        assert len(sub)==extraction['source_units'][name]
        np.testing.assert_allclose(np.asarray(raw_cache[(external.reference==name).to_numpy()].sum(axis=0)).ravel(),profiles[name].loc[axis],rtol=2e-7,atol=1e-5)
        selected=np.sort(rng.choice(sub.index.to_numpy(),5,replace=False))
        rows=external.source_row.iloc[selected].to_numpy(dtype=int)
        if name=='K562':
            path=ROOT/'reference/K562_gwps_raw_bulk_01.h5ad'
            assert digest(path)==extraction['sources'][name]['sha256']
            with h5py.File(path,'r') as f:
                ds=f['var/gene_name'];native_genes=f[ds.attrs['categories']].asstr()[:][ds[:]]
                raw=f['X'][rows,:].astype(np.float64)*f['obs/num_cells_filtered'][rows][:,None]
        elif name in ['HCT116','HEK293T']:
            path=ROOT/f'reference/orion_pilot/{name}_Batch1.parquet'
            assert digest(path)==extraction['sources'][name]['sha256']
            assert digest(ROOT/'reference/orion_pilot/gene_metadata.parquet')==extraction['sources'][name]['metadata_sha256']
            gene_meta=pq.read_table(ROOT/'reference/orion_pilot/gene_metadata.parquet').to_pandas()
            assert set(axis).issubset(set(gene_meta.gene_name))
            lookup=gene_meta.set_index('gene_token_id').gene_name
            source=pq.read_table(path,columns=['gene_token_id','gene_expression'])
            independent=[]
            for row in rows:
                tokens=source['gene_token_id'][int(row)].as_py();values=source['gene_expression'][int(row)].as_py()
                series=pd.Series(values,index=lookup.loc[tokens].to_numpy()).groupby(level=0).sum()
                # Every axis symbol is measured; sparse omissions are genuine zero counts.
                independent.append(series.reindex(axis,fill_value=0).to_numpy())
            raw=np.array(independent);native_genes=axis
        elif name=='H1':
            path=ROOT/'reference/h1_response/selected_counts.npz'
            assert digest(path)==extraction['sources'][name]['counts_sha256']
            native_genes=pd.read_csv(path.parent/'gene_names.csv').gene_name.to_numpy()
            source_meta=pd.read_csv(path.parent/'selected_cells.csv')
            pos=source_meta.set_index('source_row').index.get_indexer(rows)
            raw=sparse.load_npz(path)[pos].toarray()
        else:
            reader=CachedRangeReader(PROFILE_DIR/'ranges_CD4',extraction['sources'][name]['identity'])
            with h5py.File(reader,'r') as f:
                native_genes=column(f['var/gene_name']);x=f['X'];pointer=x['indptr'][:]
                raw=[]
                for row in rows:
                    lo,hi=int(pointer[row]),int(pointer[row+1])
                    raw.append(np.bincount(x['indices'][lo:hi],weights=x['data'][lo:hi],minlength=len(native_genes)))
                raw=np.array(raw)
            reader.close()
        matched=pd.DataFrame(raw,columns=native_genes).T.groupby(level=0,sort=False).sum().T.reindex(columns=axis)
        assert not matched.isna().any().any()
        np.testing.assert_allclose(matched,raw_cache[selected].toarray(),rtol=1e-7,atol=1e-7)
        verify_projected(matched.to_numpy(),external.iloc[selected],model)
        source_pairs.extend(matched.to_numpy()[:2])
    raw_pooled=pd.read_csv(OUT/'pooled_raw_counts.csv')
    np.testing.assert_array_equal(raw_pooled.gene_name,axis)
    raw=raw_pooled[pooled.profile_name].to_numpy().T
    verify_projected(raw,pooled,model)
    for i,row in pooled.iterrows():
        values=[]
        for j,c in enumerate(CONTEXTS):
            rho=spearmanr(raw[i],raw[j]).statistic
            assert abs(rho-row['pooled_profile_spearman_'+c])<2e-6
            values.append(rho)
        assert CONTEXTS[int(np.argmax(values))]==row.nearest_pooled_context
    for i in range(len(source_pairs)-1):
        a,b=source_pairs[i],source_pairs[i+1]
        z=rankdata(a);w=rankdata(b);z-=z.mean();w-=w.mean()
        assert abs((z@w)/(np.linalg.norm(z)*np.linalg.norm(w))-spearmanr(a,b).statistic)<1e-12
    result={'status':'passed','competition_ids_checked':55200,'external_ids_and_source_rows_checked':4393,
            'competition_raw_cells_spot_checked':30,'external_original_units_spot_checked':25,
            'pooled_profiles_independently_checked':11,'pooled_spearman_checks':33,
            'network_bytes_downloaded':0,
            'checks':'Hashes, common gene axis, all IDs, full summary fractions, original-source counts, independently reranked frozen coordinates/assignments, mean correlations, distances and pooled-profile Spearman'}
    (OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':run()
