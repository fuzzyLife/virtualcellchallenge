"""Independently validate V5 public source counts, effect arrays and all generated cell depths."""
from pathlib import Path
import argparse,json,gc
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import build_submission as build
from extract_external_mapping_controls import CachedRangeReader
from prepare_reference_controls import column,index
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v5_atlas_calibrated';A=OUT/'analysis';DATA=ROOT/'reference/v5_h1'


def validate_source():
    provenance=json.loads((DATA/'extraction.json').read_text());meta=pd.read_csv(DATA/'selected_cells.csv');x=sparse.load_npz(DATA/'selected_counts.npz')
    genes=pd.read_csv(DATA/'gene_names.csv').gene_name.to_numpy();assert meta.source_cell_id.is_unique and x.shape==(len(meta),len(genes))
    for name,key in [('selected_counts.npz','counts_sha256'),('selected_cells.csv','selected_cells_sha256'),('gene_names.csv','gene_names_sha256')]:assert build.sha256(DATA/name)==provenance[key]
    assert build.sha256(ROOT/'extract_v5_h1.py')==provenance['script_sha256']
    controls=meta[meta.target_gene=='non-targeting'];assert len(controls)==1536 and set(controls.public_split)=={'train'}
    training=x[controls.index.to_numpy()];verified=0
    for split,record in zip(['validation','test'],provenance['additional_splits']):
        folder=DATA/split;sm=pd.read_csv(folder/'selected_cells.csv');sx=sparse.load_npz(folder/'selected_counts.npz')
        mapping=sm.reset_index().set_index('source_cell_id')['index'].loc[controls.source_cell_id].to_numpy()
        assert (sx[mapping]-training).nnz==0
        reader=CachedRangeReader(ROOT/'reference/upstream_review'/('ranges_H1_'+split),record['source'],budget=1024**3)
        with h5py.File(reader,'r') as f:
            assert np.array_equal(index(f['var']),genes);labels=column(f['obs/target_gene']);identifiers=index(f['obs']);indptr=f['X/indptr'][:]
            perturb=sm[sm.target_gene!='non-targeting'].groupby('target_gene',sort=True).head(2)
            chosen=pd.concat([perturb,sm[sm.target_gene=='non-targeting'].iloc[np.linspace(0,1535,8).astype(int)]])
            for r in chosen.itertuples():
                source_row=int(r.source_row);assert identifiers[source_row]==r.source_cell_id and labels[source_row]==r.target_gene
                lo,hi=int(indptr[source_row]),int(indptr[source_row+1])
                original=np.bincount(f['X/indices'][lo:hi],weights=f['X/data'][lo:hi],minlength=len(genes))
                np.testing.assert_array_equal(original,sx.getrow(r.Index).toarray().ravel());verified+=1
        reader.close();del sx;gc.collect()
    print('Verified new original source rows:',verified,'and all duplicated selected control counts',flush=True)
    return meta,x,genes,provenance,verified


def validate_contrasts(meta,x,genes):
    with np.load(A/'h1_contrasts.npz') as z:stored={k:z[k].copy() for k in z.files}
    targets=stored['targets'];library=np.asarray(x.sum(axis=1)).ravel().astype(float);norm=x.astype(float).multiply((1e6/library)[:,None]).tocsr()
    batches=sorted(meta.batch.unique());part=meta.batch.map({b:i%2 for i,b in enumerate(batches)}).to_numpy();controls=meta.target_gene.to_numpy()=='non-targeting';batch=meta.batch.to_numpy()
    def summarize(rows):
        xx=x[rows].astype(float);ll=library[rows];n=len(rows);p=np.asarray(xx.sum(axis=0)).ravel()/ll.sum()
        residual=np.asarray(xx.power(2).sum(axis=0)).ravel()-2*p*np.asarray(xx.T.dot(ll)).ravel()+p*p*np.dot(ll,ll)
        variance=np.maximum(residual,0)/n/max(n-1,1)/(ll.mean()**2)
        return p,variance
    bases={}
    for b in batches:
        rows=np.flatnonzero(controls&(batch==b));sub=norm[rows];cm=np.asarray(sub.mean(axis=0)).ravel()
        cv=np.maximum(np.asarray(sub.power(2).mean(axis=0)).ravel()-cm**2,0)/max(len(rows)-1,1)
        bp,bv=summarize(rows);bases[b]=(cm,cv,bp,bv)
    checked=0
    for pi,partition in enumerate([0,1,None]):
        for ti,t in enumerate(targets):
            selected=meta.target_gene.to_numpy()==t
            if partition is not None:selected&=part==partition
            rows=np.flatnonzero(selected);n=len(rows);sub=norm[rows];mean=np.asarray(sub.mean(axis=0)).ravel()
            var=np.maximum(np.asarray(sub.power(2).mean(axis=0)).ravel()-mean**2,0)/max(n-1,1)
            cm=np.zeros(len(genes));cv=np.zeros(len(genes));basep=np.zeros(len(genes));basevar=np.zeros(len(genes))
            for b,count in meta.iloc[rows].batch.value_counts().items():
                m,v,bp,bv=bases[b];cm+=count/n*m;cv+=(count/n)**2*v
                fraction=library[rows[batch[rows]==b]].sum()/library[rows].sum();basep+=fraction*bp;basevar+=fraction**2*bv
            values=np.clip(np.log2((mean+1)/(cm+1)),-3,3);errors=np.sqrt(var/(mean+1)**2+cv/(cm+1)**2)/np.log(2)
            bp,bv=summarize(rows);delta=np.clip(np.log1p(50000*bp)-np.log1p(50000*basep),-3,3)
            bse=50000*np.sqrt(bv/(1+50000*bp)**2+basevar/(1+50000*basep)**2)
            for key,value in [('values',values),('errors',errors),('baselines',cm),('bulk_values',delta),('bulk_errors',bse)]:
                np.testing.assert_allclose(stored[str(pi)+'_'+key][ti],value,atol=2e-8,rtol=1e-6)
            assert stored[str(pi)+'_n'][ti]==n;checked+=5*len(genes)
    batch_record=json.loads((A/'h1_contrasts.json').read_text())['globally_disjoint_batch_controls']
    assert not set(batch_record['left_control_ids'])&set(batch_record['right_control_ids'])
    print('Recalculated public H1 contrasts and pooled uncertainty:',checked,'values',flush=True)
    return stored,checked


def validate_model(stored,genes):
    manifest,official,targets=build.read_spec(ROOT);record=json.loads((A/'reference_model.json').read_text())
    assert record['script_sha256']==build.sha256(ROOT/'fit_v5_submission.py')
    assert record['reused_v4_fitter_sha256']==build.sha256(ROOT/'fit_lineage_submission.py')
    assert record['model_arrays_sha256']==build.sha256(A/'reference_model_arrays.npz')
    assert record['h1_components_sha256']==build.sha256(A/'h1_components.npz')
    for name,key in [('h1_contrasts.npz','h1_contrasts_sha256'),('h1_contrasts.json','h1_contrasts_metadata_sha256'),('promoter_selection.json','promoter_selection_sha256')]:assert build.sha256(A/name)==record[key]
    assert build.sha256(DATA/'extraction.json')==record['h1_source_extraction_sha256']
    parent=ROOT/'experiments/v4_lineage_response/analysis';old=json.loads((parent/'reference_model.json').read_text())
    assert build.sha256(parent/'reference_model_arrays.npz')==old['model_arrays_sha256']==record['v4_parent_arrays_sha256']
    assert build.sha256(parent/'source_components.npz')==old['source_components_sha256']==record['v4_source_components_sha256']
    for path,digest in old['source_data_sha256'].items():assert build.sha256(ROOT/path)==digest,path
    with np.load(parent/'source_components.npz') as z:oldc={k:z[k].copy() for k in z.files}
    with np.load(parent/'reference_model_arrays.npz') as z:old_effect=z['effects'].copy()
    with np.load(A/'reference_model_arrays.npz') as z:effect=z['effects'].copy();bulk=z['bulk_delta'].copy()
    with np.load(A/'h1_components.npz') as z:hc={k:z[k].copy() for k in z.files}
    ht=stored['targets'];hg=stored['genes'];gp={g:i for i,g in enumerate(official)};si=np.array([i for i,g in enumerate(hg) if g in gp]);di=np.array([gp[hg[i]] for i in si])
    v=stored['2_values'];se=stored['2_errors'];base=stored['2_baselines'];bval=stored['2_bulk_values'];bse=stored['2_bulk_errors']
    own=np.array([int(np.flatnonzero(hg==t)[0]) for t in ht]);ownv=v[np.arange(len(ht)),own];ownse=se[np.arange(len(ht)),own]
    accepted=(stored['2_n']>=40)&(ownv<=np.log2(.7))&(base[np.arange(len(ht)),own]>5)&(ownv<-1.96*ownse)
    valid=(base>5)&(se>0)&(hg[None,:]!=ht[:,None])&accepted[:,None]
    settings=record['h1_settings'];transfer=record['cross_lineage_settings'];checks=0
    for i,t in enumerate(targets):
        train=ht!=t;idx=np.flatnonzero(ht==t);ok=bool(len(idx) and accepted[idx[0]])
        gain=settings['prior_gain'] if ok else transfer['prior_gain']
        for prefix,values,errors,mask in [('h1',v,se,valid),('h1_bulk',bval,bse,valid&(bse>0))]:
            weight=mask[train]/(.04+errors[train]**2);den=weight.sum(axis=0)
            expected_prior=gain*np.clip(np.divide((weight*values[train]).sum(axis=0),den,out=np.zeros(len(hg)),where=den>0),-.5,.5)
            expected=np.zeros(len(hg))
            if ok:
                row=idx[0];coefficient=settings['specific_gain'] if prefix=='h1' else record['bulk_gain']
                expected=np.clip(coefficient*np.sign(values[row])*np.maximum(np.abs(values[row])-errors[row],0),-3,3);expected[~valid[row]]=0
            np.testing.assert_allclose(hc[prefix][i,di],expected[si],atol=2e-7,rtol=2e-6)
            np.testing.assert_allclose(hc[prefix+'_prior'][i,di],expected_prior[si],atol=2e-7,rtol=2e-6);checks+=2*len(di)
        assert bool(hc['h1_available'][i])==ok
        expected_mask=np.zeros(len(official),dtype=bool)
        if ok:expected_mask[di]=valid[idx[0],si]
        np.testing.assert_array_equal(hc['h1_mask'][i],expected_mask)
    # A inherits the frozen V4 effect except explicitly recorded promoter changes.
    changed=pd.read_csv(A/'promoter_effect_changes.csv');expectedA=old_effect[0].copy();ti={t:i for i,t in enumerate(targets)}
    selection=json.loads((A/'promoter_selection.json').read_text())
    assert selection['script_sha256']==build.sha256(ROOT/'calibrate_v5_promoters.py')
    for path,digest in selection['source_sha256'].items():assert build.sha256(ROOT/path)==digest,path
    pairs=pd.read_csv(ROOT/'analysis/upstream_review/promoter_neighbor_pairs.csv')
    pm=pd.read_csv(ROOT/'reference/cd4_probe/merged.csv')
    with np.load(ROOT/'reference/cd4_probe/merged.npz') as z:pv=z['log_fc'].copy();pq=z['adj_p_value'].copy();pn=z['genes'].copy()
    pg={g:i for i,g in enumerate(pn)}
    quality=(pm.n_cells_target>=100)&pm.ontarget_significant&~pm.distal_offtarget_flag&~pm.single_guide_estimate&~pm.low_target_gex
    expected_changes=[]
    if selection['accepted']:
        for pair in pairs.itertuples():
            if pair.neighbor not in pg:continue
            rows=np.flatnonzero(quality.to_numpy()&(pm.target_contrast_gene_name.to_numpy()==pair.target));j=pg[pair.neighbor]
            support=(pv[rows,j]<0)&(pq[rows,j]<.05)&np.isfinite(pv[rows,j])&np.isfinite(pq[rows,j])
            if pm.iloc[rows[support]].culture_condition.nunique()<2:continue
            i,j=ti[pair.target],gp[pair.neighbor]
            scale=1-np.clip(np.log(max(pair.distance,500)/500)/np.log(10),0,1)
            gg=min(1,oldc['context_mean_cpm'][0,j]/5);tg=min(1,oldc['context_mean_cpm'][0,gp[pair.target]]/5)
            value=min(float(expectedA[i,j]),selection['fitted_close_neighbor_log2fc']*scale*gg*tg)
            if np.float32(value)!=expectedA[i,j]:expected_changes.append((pair.target,pair.neighbor))
            expectedA[i,j]=value
    assert set(expected_changes)==set(zip(changed.target_gene,changed.neighbor))
    np.testing.assert_array_equal(effect[0],expectedA)
    for ci,c in enumerate('ABC'):
        obj=ad.read_h5ad(ROOT/f'context_{c}.h5ad');x=obj.X.tocsr();p=np.asarray(x.sum(axis=0)).ravel().astype(float);p/=p.sum();basebulk=50000*p
        if ci==0:expected_bulk=np.log1p(basebulk[None,:]*np.exp2(effect[0]))-np.log1p(basebulk[None,:])
        else:
            source=transfer['cd4_coefficient']*oldc['cd4_specific'][ci]+transfer['k562_coefficient']*oldc['k562']
            gg=np.clip(oldc['context_mean_cpm'][ci]/5,0,1);tg=np.clip(oldc['context_mean_cpm'][ci,[gp[t] for t in targets]]/5,0,1)
            expected=np.clip((np.where(hc['h1_mask'],hc['h1'],source)+hc['h1_prior'])*gg[None,:]*tg[:,None],-3,3)
            for i,t in enumerate(targets):expected[i,gp[t]]=0
            np.testing.assert_allclose(effect[ci],expected,atol=2e-7,rtol=3e-6)
            translated=np.log1p(basebulk[None,:]*np.exp2(np.clip(source,-3,3)))-np.log1p(basebulk[None,:])
            expected_bulk=(np.where(hc['h1_mask'],hc['h1_bulk'],record['bulk_gain']*translated)+hc['h1_bulk_prior'])*gg[None,:]*tg[:,None]
        for i,t in enumerate(targets):expected_bulk[i,gp[t]]=0
        np.testing.assert_allclose(bulk[ci],expected_bulk,atol=3e-7,rtol=3e-6);obj.file.close();del obj,x;gc.collect()
    for filename in ['h1_lineage_selection.json','bulk_selection.json','promoter_selection.json']:
        selection=json.loads((A/filename).read_text());seen=[]
        for fold in selection['nested_folds']:
            assert not set(fold['train_targets'])&set(fold['test_targets']);seen+=fold['test_targets']
        assert len(seen)==len(set(seen))
    print('Independently assembled all mean and pooled V5 effects',flush=True)
    return record,effect,bulk,checks


def validate_package(record,effects,bulk):
    from build_v5_submission import emit
    from evaluate_atlas_transfer import load_atlas_model
    manifest,genes,targets=build.read_spec(ROOT);rows=pd.read_csv(A/'prediction_control_rows.csv');stats=pd.read_csv(A/'prediction_group_statistics.csv')
    assert len(rows)==360000 and len(stats)==900 and (rows.groupby(['context','target_gene']).size()==400).all()
    params=json.loads((A/'run_parameters.json').read_text());assert params['generator_sha256']==build.sha256(ROOT/'build_v5_submission.py')
    assert params['source_model_sha256']==build.sha256(A/'reference_model.json')
    package=json.loads((A/'submission_ready.json').read_text());assert package['official_prep_passed'] and package['bytes']==(OUT/'prediction.vcc').stat().st_size
    assert package['sha256']==build.sha256(OUT/'prediction.vcc')
    predicted=ad.read_h5ad(OUT/'prediction.h5ad',backed='r');assert predicted.shape==(360000,18533) and np.array_equal(predicted.var_names,genes)
    counts_checked=0;replayed=0;model,_=load_atlas_model()
    for ci,c in enumerate('ABC'):
        obj=ad.read_h5ad(ROOT/f'context_{c}.h5ad');x=obj.X.tocsr()
        for start in range(0,x.nnz,4000000):
            values=x.data[start:start+4000000]
            assert np.isfinite(values).all() and (values>=0).all() and (values==np.floor(values)).all()
        # Exact replay uses the same canonical integer counts as generation, avoiding float32 column-sum rounding.
        x=x.astype(np.int32);x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
        library=np.asarray(x.sum(axis=1,dtype=np.int64)).ravel();guide=obj.obs.ntc_id.astype(str).to_numpy()
        source=rows[rows.context==c];indices=source.source_control_row.to_numpy();assert ((indices>=0)&(indices<len(library))).all()
        np.testing.assert_array_equal(source.original_total_counts,library[indices]);np.testing.assert_array_equal(source.ntc_id,guide[indices]);np.testing.assert_array_equal(source.source_control_cell_id,obj.obs_names.to_numpy()[indices])
        for target,group in source.groupby('target_gene'):
            assert group.source_control_row.nunique()==400
            ti=list(targets).index(target);rng=np.random.default_rng(np.random.SeedSequence([params['seed'],ci,ti]))
            expected_rows=build.guide_stratified_sample(guide,library,400,rng)
            np.testing.assert_array_equal(group.source_control_row.to_numpy(),expected_rows)
            guide_counts=group.ntc_id.value_counts();assert len(guide_counts)==46 and guide_counts.max()-guide_counts.min()==1
        # Exact generator replay uses its documented 512-row accumulation order; source moments were independently checked above.
        cpm_mean=np.zeros(len(genes),dtype=float)
        for first in range(0,x.shape[0],512):
            sub=x[first:first+512].astype(float)
            cpm_mean+=np.asarray(sub.multiply((1e6/library[first:first+512])[:,None]).sum(axis=0)).ravel()
        cpm_mean/=x.shape[0];mean=cpm_mean/cpm_mean.sum();pooled=np.asarray(x.sum(axis=0)).ravel().astype(float);pooled/=pooled.sum()
        for start in range(ci*120000,(ci+1)*120000,1000):
            block=predicted.X[start:start+1000].tocsr();assert np.issubdtype(block.dtype,np.integer) and (block.data>0).all()
            expected=source.iloc[start-ci*120000:start-ci*120000+1000].original_total_counts.to_numpy()
            np.testing.assert_array_equal(np.asarray(block.sum(axis=1)).ravel(),expected);counts_checked+=block.nnz
        for target in ['ACLY','ALG12','MED25','VIM']:
            ti=list(targets).index(target);group=source[source.target_gene==target];idx=group.source_control_row.to_numpy()
            p=mean*np.exp2(np.clip(effects[ci,ti],-3,3));b=np.expm1(np.maximum(np.log1p(50000*pooled)+np.clip(bulk[ci,ti],-3,3),0));own=list(genes).index(target)
            p[own]*=1-params['knockdown_fraction_assumption'];b[own]*=1-params['knockdown_fraction_assumption'];p/=p.sum();b/=b.sum()
            block,_=emit(model,x,library,idx,p,b,params['seed']+ci*100000+ti)
            actual=predicted.X[(ci*300+ti)*400:(ci*300+ti+1)*400].tocsr();assert (actual-block).nnz==0;replayed+=1
        obj.file.close();del obj,x;gc.collect()
    predicted.file.close()
    assert counts_checked==package['prediction_validation']['nnz']
    assert stats.requested_mean_probability_L1.max()<=.03 and stats.projected_pooled_probability_L1.max()<=.002
    assert stats.bulk_feasibility_projection_L1.max()<=.03
    print('All cells/counts/depths checked; fully replayed groups:',replayed,flush=True)
    return {'all_original_control_records':len(rows),'all_integer_count_entries':counts_checked,'exact_generated_groups_replayed':replayed,'package_sha256':package['sha256'],
            'maximum_mean_probability_L1':float(stats.requested_mean_probability_L1.max()),'maximum_projected_pooled_probability_L1':float(stats.projected_pooled_probability_L1.max())}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--model-only',action='store_true');args=parser.parse_args()
    meta,x,genes,provenance,original_checks=validate_source();contrasts,contrast_checks=validate_contrasts(meta,x,genes);del x;gc.collect()
    record,effects,bulk,effect_checks=validate_model(contrasts,genes)
    evidence={'passed':True,'original_public_H1_rows_independently_checked':original_checks,'repeated_controls_exact_count_checked_per_split':1536,
              'source_contrast_values_recalculated':contrast_checks,'H1_direct_and_prior_values_reassembled':effect_checks,
              'all_mean_and_bulk_effects_reassembled':int(effects.size+bulk.size),'script_sha256':build.sha256(Path(__file__)),
              'source_model_sha256':build.sha256(A/'reference_model.json'),'scope':'Public cached-source provenance, unique controls, original expression spot checks, complete source moments/effects, target/batch separation; no hidden-reference or leaderboard validation.'}
    build.write_json(A/'model_validation.json',evidence)
    if not args.model_only:
        evidence['generated_package_checks']=validate_package(record,effects,bulk)
        build.write_json(A/'independent_validation.json',evidence)
    print(json.dumps(evidence,indent=2),flush=True)


if __name__=='__main__':main()
