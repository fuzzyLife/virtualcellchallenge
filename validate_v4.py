"""Independently validate V4 source calibration, effect assembly and generated package."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import anndata as ad
import h5py
from scipy import sparse
import build_submission as build

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'experiments/v4_lineage_response';A=OUT/'analysis'

def source_checks(model,genes,targets,components):
    # Recalculate K562 response means from original rows, independently of the fitter.
    with h5py.File(ROOT/'reference/K562_gwps_raw_bulk_01.h5ad','r') as f:
        labels=f['obs/gene_transcript'].asstr()[:];number=f['obs/num_cells_filtered'][:].astype(float)
        reduction=f['obs/fold_expr'][:].astype(float);core=f['obs/core_control'][:].astype(bool)
        column=f['var/gene_name'];symbols=f[column.attrs['categories']].asstr()[:][column[:]]
        source_target=np.array([x.split('_',2)[1] if x.split('_',1)[0].isdigit() and len(x.split('_',2))==3 else '' for x in labels])
        shared=sorted(set(symbols)&set(genes));position={g:i for i,g in enumerate(shared)}
        idx=np.array([i for i,g in enumerate(symbols) if g in position]);dest=np.array([position[symbols[i]] for i in idx])
        projection=sparse.csr_matrix((np.ones(len(idx)),(idx,dest)),shape=(len(symbols),len(shared)))
        control=np.flatnonzero(core&(source_target=='non-targeting')&np.isfinite(number)&(number>=30))
        cm=np.asarray(projection.T.dot(f['X'][control,:].astype(float).T).T)
        basal=np.average(cm,axis=0,weights=number[control]);bcpm=basal*1e6/basal.sum();cn=number[control].sum()
        effective=np.isfinite(number)&(number>=30)&np.isfinite(reduction)&(reduction<=.5)&(source_target!='non-targeting')
        gp=np.array([list(genes).index(g) for g in shared]);checks=0
        for ti,target in enumerate(targets):
            rows=np.flatnonzero(effective&(source_target==target))
            expected=np.zeros(len(shared))
            if len(rows):
                values=np.asarray(projection.T.dot(f['X'][rows,:].astype(float).T).T)
                mean=np.average(values,axis=0,weights=number[rows]);cpm=mean*1e6/mean.sum();n=number[rows].sum()
                raw=np.log2((cpm+1)/(bcpm+1));se=np.sqrt(1/np.maximum(n*mean,1)+1/np.maximum(cn*basal,1))/np.log(2)
                expected=.3*np.clip(np.sign(raw)*np.maximum(np.abs(raw)-1.96*se,0),-3,3)*n/(n+50)
            np.testing.assert_allclose(components['k562'][ti,gp],expected,atol=3e-8,rtol=3e-6)
            checks+=len(shared)
    print('Original K562 means independently checked for all 300 targets',flush=True)
    # Recalculate mapped CD4 target effects on original merged rows for all measured genes.
    meta=pd.read_csv(ROOT/'reference/cd4_probe/merged.csv')
    with np.load(ROOT/'reference/cd4_probe/merged.npz',allow_pickle=False) as z:
        source=z['genes'];raw=np.clip(np.nan_to_num(z['log_fc'].astype(float),nan=0,posinf=0,neginf=0),-3,3)
        se=np.nan_to_num(z['lfcSE'].astype(float),nan=1e6,posinf=1e6,neginf=1e6)
        valid=np.isfinite(z['log_fc'])&np.isfinite(z['lfcSE'])&(z['lfcSE']>0)&(source[None,:]!=meta.target_contrast_gene_name.to_numpy()[:,None])
    good=((meta.n_cells_target>=100)&meta.ontarget_significant.astype(bool)&~meta.distal_offtarget_flag.astype(bool)&~meta.neighboring_gene_KD.astype(bool)&~meta.single_guide_estimate.astype(bool)&~meta.low_target_gex.astype(bool)).to_numpy()
    positions={g:i for i,g in enumerate(genes)};si=np.array([i for i,g in enumerate(source) if g in positions]);di=np.array([positions[source[i]] for i in si])
    condition_names=['Rest','Stim8hr','Stim48hr'];setting=model['cd4_settings'];checks_cd4=0
    for ti,target in enumerate(targets):
        specific=np.zeros((3,len(si)));masks=np.zeros_like(specific,dtype=bool);prior=np.zeros_like(specific)
        for ci,name in enumerate(condition_names):
            use=np.flatnonzero(good&(meta.target_contrast_gene_name.to_numpy()==target)&(meta.culture_condition.to_numpy()==name))
            if len(use):
                weight=np.sqrt(meta.n_cells_target.to_numpy()[use,None])*valid[use][:,si]
                total=weight.sum(axis=0)
                estimate=setting['specific_gain']*raw[use][:,si]*.04/(.04+se[use][:,si]**2)
                specific[ci]=np.divide((weight*estimate).sum(axis=0),total,out=np.zeros(len(si)),where=total>0);masks[ci]=total>0
            train=np.flatnonzero(good&(meta.target_contrast_gene_name.to_numpy()!=target)&(meta.culture_condition.to_numpy()==name))
            weight=valid[train][:,si]/(se[train][:,si]**2+.04);total=weight.sum(axis=0)
            prior[ci]=setting['prior_gain']*np.clip(np.divide((weight*raw[train][:,si]).sum(axis=0),total,out=np.zeros(len(si)),where=total>0),-.5,.5)
        for ci,context in enumerate(['A','B','C']):
            weights=np.array([model['context_mixture_weights'][context]['cd4_conditions']['CD4_'+name] for name in condition_names]);weights=.5*weights+.5/3
            weight=weights[:,None]*masks;den=weight.sum(axis=0)
            expected=np.divide((weight*specific).sum(axis=0),den,out=np.zeros(len(si)),where=den>0)
            np.testing.assert_allclose(components['cd4_specific'][ci,ti,di],expected,atol=1e-7,rtol=3e-6)
            expected=np.clip(expected+np.einsum('c,cg->g',weights,prior),-3,3)
            np.testing.assert_allclose(components['cd4'][ci,ti,di],expected,atol=1e-7,rtol=3e-6)
            np.testing.assert_array_equal(components['cd4_mask'][ci,ti,di],masks.any(axis=0));checks_cd4+=len(si)
    print('Original CD4 effects independently checked for every target/context',flush=True)
    # Recalculate H1 full-source effects and leave-own-target-out shared means.
    path=ROOT/'reference/h1_response';x=sparse.load_npz(path/'selected_counts.npz').astype(float)
    meta=pd.read_csv(path/'selected_cells.csv');source=pd.read_csv(path/'gene_names.csv').gene_name.to_numpy()
    libraries=np.asarray(x.sum(axis=1)).ravel();cpm=x.multiply((1e6/libraries)[:,None]).tocsr()
    htargets=sorted(set(meta.target_gene)-{'non-targeting'});bs={}
    for batch in sorted(meta.batch.unique()):
        rows=np.flatnonzero((meta.batch.to_numpy()==batch)&(meta.target_gene.to_numpy()=='non-targeting'))
        sub=cpm[rows];mean=np.asarray(sub.mean(axis=0)).ravel();variance=np.maximum(np.asarray(sub.power(2).mean(axis=0)).ravel()-mean**2,0)/max(len(rows)-1,1)
        bs[batch]=(mean,variance)
    values=[];errors=[];supported=[];nvalues=[]
    for target in htargets:
        rows=np.flatnonzero(meta.target_gene.to_numpy()==target);sub=cpm[rows];n=len(rows);mean=np.asarray(sub.mean(axis=0)).ravel()
        vm=np.maximum(np.asarray(sub.power(2).mean(axis=0)).ravel()-mean**2,0)/max(n-1,1)
        base=np.zeros(len(source));bv=np.zeros(len(source))
        for batch,num in meta.iloc[rows].batch.value_counts().items():
            bm,bvar=bs[batch];base+=num/n*bm;bv+=(num/n)**2*bvar
        lfc=np.clip(np.log2((mean+1)/(base+1)),-3,3);se=np.sqrt(vm/(mean+1)**2+bv/(base+1)**2)/np.log(2)
        values.append(lfc);errors.append(se);supported.append((base>5)&(se>0)&(source!=target));nvalues.append(n)
    values=np.array(values);errors=np.array(errors);supported=np.array(supported);si=np.array([i for i,g in enumerate(source) if g in positions]);di=np.array([positions[source[i]] for i in si])
    hs=model['h1_settings'];transfer=model['cross_lineage_settings'];checks_h1=0
    for ti,target in enumerate(targets):
        own=np.array(htargets)==target;train=~own;weight=supported[train]/(errors[train]**2+.04);den=weight.sum(axis=0)
        mean=np.clip(np.divide((weight*values[train]).sum(axis=0),den,out=np.zeros(len(source)),where=den>0),-.5,.5)
        direct=np.zeros(len(source));accepted=False;mask=np.zeros(len(source),bool)
        if own.any():
            j=np.flatnonzero(own)[0];gi=np.flatnonzero(source==target)
            accepted=nvalues[j]>=40 and len(gi)>0 and values[j,gi[0]]<=np.log2(.7)
            if accepted:
                mask=supported[j];direct=hs['specific_gain']*np.sign(values[j])*np.maximum(np.abs(values[j])-errors[j],0);direct[~mask]=0
        background=(hs['prior_gain'] if accepted else transfer['prior_gain'])*mean
        background[~supported.any(axis=0)]=0
        np.testing.assert_allclose(components['h1'][ti,di],direct[si],rtol=3e-6,atol=1e-7)
        np.testing.assert_allclose(components['h1_prior'][ti,di],background[si],rtol=3e-6,atol=1e-7)
        np.testing.assert_array_equal(components['h1_mask'][ti,di],mask[si]);assert bool(components['h1_available'][ti])==bool(accepted)
        checks_h1+=2*len(si)
    print('Original H1 direct and held-own-target-out priors independently checked',flush=True)
    return {'k562_effect_values':checks,'cd4_context_effect_values':checks_cd4,'h1_direct_and_prior_values':checks_h1}


def model_check():
    manifest,genes,targets=build.read_spec(ROOT);model=json.loads((A/'reference_model.json').read_text())
    assert build.sha256(ROOT/'fit_lineage_submission.py')==model['script_sha256']
    assert build.sha256(ROOT/'fit_context_ensemble.py')==model['reuse_fitter_sha256']
    assert build.sha256(A/'reference_model_arrays.npz')==model['model_arrays_sha256']
    assert build.sha256(A/'source_components.npz')==model['source_components_sha256']
    for path,value in model['source_data_sha256'].items():assert build.sha256(ROOT/path)==value,path
    for record in model['context_original_inputs']:assert build.sha256(ROOT/record['file'])==record['sha256']
    with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as z:
        effects=z['effects'].copy();np.testing.assert_array_equal(z['genes'],genes);np.testing.assert_array_equal(z['targets'],targets)
        np.testing.assert_array_equal(z['cell_lines'],['A','B','C']);assert effects.shape==(3,300,18533)
        np.testing.assert_array_equal(z['available'],np.any(effects!=0,axis=2))
    with np.load(A/'source_components.npz',allow_pickle=False) as z:components={key:z[key].copy() for key in z.files}
    assert np.isfinite(effects).all() and np.abs(effects).max()<=3
    positions={g:i for i,g in enumerate(genes)};own=np.array([positions[t] for t in targets])
    transfer=model['cross_lineage_settings'];source_evidence=source_checks(model,genes,targets,components)
    for ci,context in enumerate(['A','B','C']):
        assert build.sha256(ROOT/f'analysis/gene_statistics_{context}.csv')==model['context_statistics_sha256'][context]
        stats_frame=pd.read_csv(ROOT/f'analysis/gene_statistics_{context}.csv');np.testing.assert_array_equal(stats_frame.gene_name,genes)
        means=stats_frame.mean_cpm.to_numpy();np.testing.assert_array_equal(means,components['context_mean_cpm'][ci])
        if ci==0:expected=np.where(components['cd4_mask'][ci],components['cd4'][ci],np.where(components['k562_mask'],components['k562'],0))
        else:
            translated=transfer['cd4_coefficient']*components['cd4_specific'][ci]+transfer['k562_coefficient']*components['k562']
            expected=np.where(components['h1_mask'],components['h1'],translated)+components['h1_prior']
        expected=np.clip(expected*np.clip(means/5,0,1)[None,:]*np.clip(means[own]/5,0,1)[:,None],-3,3)
        expected[np.arange(len(targets)),own]=0
        np.testing.assert_allclose(effects[ci],expected,atol=1e-7,rtol=3e-6)
    all_folds=[]
    for key in ['cd4_selection','h1_selection','cross_lineage_selection']:
        selection=model[key];test=[]
        for fold in selection['nested_folds']:
            assert not set(fold['train_targets'])&set(fold['test_targets']);test.extend(fold['test_targets'])
        assert len(test)==len(set(test))==selection['targets']
        folds=selection['nested_folds'];n=sum(x['contrasts'] for x in folds)
        for column in ['loss','zero_loss','selection_score','direction_fidelity_proxy','jaccard_proxy','mean_predicted_calls','mean_reference_calls']:
            expected=sum(x[column]*x['contrasts'] for x in folds)/n
            np.testing.assert_allclose(selection['nested_procedure'][column],expected,atol=1e-12)
        np.testing.assert_allclose(selection['nested_procedure']['relative_zero_loss'],selection['nested_procedure']['loss']/selection['nested_procedure']['zero_loss'],atol=1e-12)
        all_folds.append({'source':key,'targets':len(test),'folds':len(folds),'held_target_overlap':False})
    records=pd.read_csv(A/'target_response_sources.csv');assert len(records)==900
    for ci,context in enumerate(['A','B','C']):
        part=records[records.context==context];np.testing.assert_array_equal(part.target_gene,targets)
        np.testing.assert_array_equal(part.nonzero_effect_genes,(effects[ci]!=0).sum(axis=1))
    evidence={'passed':True,'source_data_checksums_passed':True,'all_16679700_effect_values_assembled_independently':True,
       'all_source_effect_values_recomputed':source_evidence,'nested_target_fold_checks':all_folds,
       'fitter_sha256':model['script_sha256'],'validator_sha256':build.sha256(Path(__file__)),
       'scope':'Source/numerical validation and 400-cell call proxies, not a hidden-label competition score'}
    build.write_json(A/'model_validation.json',evidence);print(json.dumps(evidence,indent=2),flush=True)
    return manifest,genes,targets,effects


def prediction_check(manifest,genes,targets,effects):
    ready=json.loads((A/'submission_ready.json').read_text());assert ready['official_prep_passed'] and not ready['uploaded']
    assert build.sha256(OUT/'prediction.vcc')==ready['sha256']
    parameters=json.loads((A/'run_parameters.json').read_text());assert parameters['reference_kind']=='lineage-response'
    prediction=ad.read_h5ad(OUT/'prediction.h5ad',backed='r');assert prediction.shape==(360000,18533)
    np.testing.assert_array_equal(prediction.var_names,genes)
    assert prediction.obs_names.is_unique
    assert prediction.obs.groupby(['context','target_gene'],observed=True).size().to_dict()=={(c,t):400 for c in ['A','B','C'] for t in targets}
    source_rows=pd.read_csv(A/'prediction_control_rows.csv');assert len(source_rows)==360000
    stats_frame=pd.read_csv(A/'prediction_group_statistics.csv');assert len(stats_frame)==900
    checks=[];selected=set([0,299])
    evidence=pd.read_csv(A/'target_response_sources.csv').query("context=='C'")
    for mask in [evidence.h1_direct_response_available.to_numpy(),evidence.cd4_measured_genes.to_numpy()==0]:
        if mask.any():selected.add(int(np.flatnonzero(mask)[0]))
    energy=np.sqrt((effects.astype(float)**2).mean(axis=(0,2)))
    selected.update(np.argsort(energy)[np.rint(np.linspace(0,299,4)).astype(int)].tolist())
    for ci,context in enumerate(['A','B','C']):
        original=ad.read_h5ad(ROOT/f'context_{context}.h5ad')
        guides=original.obs.ntc_id.astype(str).to_numpy();x=original.X.tocsr().astype(np.int32)
        x.sum_duplicates();x.eliminate_zeros();x.sort_indices();libraries=np.asarray(x.sum(axis=1)).ravel()
        means=np.asarray(x.mean(axis=0)).ravel();mean_lib=float(libraries.mean())
        rows_context=source_rows[source_rows.context==context]
        for ti,target in enumerate(targets):
            entry=rows_context[rows_context.target_gene==target];rows=entry.source_control_row.to_numpy()
            assert len(rows)==len(set(rows))==400 and (rows>=0).all() and (rows<18400).all()
            np.testing.assert_array_equal(entry.prediction_cell_offset,np.arange(400))
            np.testing.assert_array_equal(entry.ntc_id,guides[rows]);np.testing.assert_array_equal(entry.original_total_counts,libraries[rows])
            counts=pd.Series(guides[rows]).value_counts();assert len(counts)==46 and counts.min()==8 and counts.max()==9
            # Reproduce every original sampling choice from the documented seed.
            rng=np.random.default_rng(np.random.SeedSequence([parameters['seed'],ci,ti]));expected_rows=build.guide_stratified_sample(guides,libraries,400,rng)
            np.testing.assert_array_equal(rows,expected_rows)
            if ti in selected:
                sample=x[rows];totals=np.asarray(sample.sum(axis=0)).ravel();start=(ci*300+ti)*400
                block=prediction.X[start:start+400].tocsr()
                assert np.issubdtype(block.dtype,np.integer) and (block.data>0).all()
                factor=np.exp2(effects[ci,ti].astype(float));factor[list(genes).index(target)]=1-parameters['knockdown_fraction_assumption']
                down=np.minimum(factor,1);up=libraries[rows].sum()/mean_lib*means*np.maximum(factor-1,0)
                expected=totals*down+up;variance=totals*down*(1-down)+up;actual=np.asarray(block.sum(axis=0)).ravel()
                assert (np.abs(actual-expected)<=8*np.sqrt(variance)+10).all(),(context,target)
                np.testing.assert_array_equal(actual[factor==1],totals[factor==1])
                # Exact selected-cell generation replay verifies origin, axis and stochastic seed.
                replay=build.simulate_counts_zero_aware(x,rows,factor,rng,means,mean_lib)
                assert (replay!=block).nnz==0,(context,target,'count replay')
                entry_stats=stats_frame[(stats_frame.context==context)&(stats_frame.target_gene==target)].iloc[0]
                assert block.nnz==entry_stats.nnz
                checks.append({'context':context,'target_gene':target,'gene_count_expectations':len(genes),'exact_generated_cell_replay':True})
        print(context,'all 300 control samples and selected count blocks checked',flush=True)
        original.file.close();del original,x
    prediction.file.close()
    from vcc.vccfile import validate_vcc,read_vcc_meta
    validate_vcc(OUT/'prediction.vcc');meta=read_vcc_meta(OUT/'prediction.vcc');assert meta['nnz']==ready['prediction_validation']['nnz']
    result={'passed':True,'official_package_validated':True,'package_sha256_matches':True,'all_900_groups_validated':True,
       'all_360000_source_control_rows_verified':True,'independent_count_expectation_groups':checks,
       'independent_gene_count_expectations':len(checks)*len(genes),'validator_sha256':build.sha256(Path(__file__)),
       'scope':'Full submission validity, source provenance, independent count expectations and deterministic replay; no hidden-label accuracy claim'}
    build.write_json(A/'independent_validation.json',result);print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--model-only',action='store_true');parser.add_argument('--prediction-only',action='store_true')
    args=parser.parse_args()
    if args.prediction_only:
        manifest,genes,targets=build.read_spec(ROOT)
        with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as z:effects=z['effects'].copy()
    else:manifest,genes,targets,effects=model_check()
    if not args.model_only:prediction_check(manifest,genes,targets,effects)
