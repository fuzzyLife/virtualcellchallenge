"""Validate the V3 ensemble and independently audit generated count expectations."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.stats import ks_2samp
import build_submission as build
from audit_prediction import ControlRankTest, normalized_counts

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v3_context_ensemble'
A=OUT/'analysis'


def validate_model():
    manifest,genes,targets=build.read_spec(ROOT)
    model=json.loads((A/'reference_model.json').read_text())
    assert build.sha256(A/'reference_model_arrays.npz')==model['model_arrays_sha256']
    assert build.sha256(ROOT/'fit_context_ensemble.py')==model['script_sha256']
    for name,checksum in model['source_data_sha256'].items():
        assert build.sha256(ROOT/'reference/cd4_probe'/name)==checksum
    with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as z:
        effects=z['effects'].copy()
        np.testing.assert_array_equal(z['genes'],genes)
        np.testing.assert_array_equal(z['targets'],targets)
        assert effects.shape==(3,300,18533) and np.isfinite(effects).all()
        assert np.abs(effects).max()<=3.00001
        np.testing.assert_array_equal(z['cell_lines'],['A','B','C'])
    assert not np.array_equal(effects[0],effects[1]) and not np.array_equal(effects[1],effects[2])
    for ti,target in enumerate(targets):
        assert (effects[:,ti,np.flatnonzero(genes==target)[0]]==0).all()
    selection=model['cd4_model_selection']
    test_targets=[]
    for fold in selection['nested_folds']:
        assert not set(fold['train_targets'])&set(fold['test_targets'])
        test_targets.extend(fold['test_targets'])
    assert len(test_targets)==len(set(test_targets))==selection['targets']
    assert selection['nested_selection_relative_zero_loss']<selection['v2_outer_cv_relative_zero_loss']
    h1=pd.read_csv(A/'h1_target_holdout_validation.csv')
    np.testing.assert_allclose(h1.loss.sum()/h1.zero_loss.sum(),model['h1_model_selection']['nested_target_holdout_relative_zero_loss'],atol=1e-12)
    weights=pd.read_csv(A/'target_source_weights.csv')
    assert len(weights)==900
    for context in ['A','B','C']:
        sub=weights[weights.context==context]
        np.testing.assert_array_equal(sub.target_gene,targets)
        np.testing.assert_array_equal(sub.target_index,np.arange(300))
        w=sub[['cd4_weight','k562_weight','h1_weight']].to_numpy()
        assert np.isfinite(w).all() and (w>=0).all()
        sums=w.sum(axis=1)
        support=sub[['cd4_available','k562_available','h1_available']].any(axis=1).to_numpy()
        np.testing.assert_allclose(sums,np.where(support,1.,0.),atol=1e-12)
    template_checks=[]
    previous=pd.read_csv(ROOT/'analysis/target_rank_correlations/per_cell_target_rank_correlations.csv')
    for ci,context in enumerate(['A','B','C']):
        template=pd.read_csv(A/f'control_template_{context}.csv')
        obj=ad.read_h5ad(ROOT/f'context_{context}.h5ad',backed='r')
        rows=template.control_row.to_numpy()
        assert len(rows)==len(set(rows))==400 and (rows>=0).all() and (rows<18400).all()
        counts=obj.X[rows].tocsr()
        actual=np.asarray(counts.sum(axis=1)).ravel()
        np.testing.assert_array_equal(actual,template.total_counts)
        np.testing.assert_array_equal(template.ntc_id,obj.obs.ntc_id.astype(str).to_numpy()[rows])
        balance=template.ntc_id.value_counts()
        assert len(balance)==46 and balance.min()==8 and balance.max()==9
        all_depth=previous[previous.context==context].total_counts.to_numpy()
        template_checks.append({'context':context,'cells':400,'guides':46,'min_cells_per_guide':8,'max_cells_per_guide':9,
            'depth_ks_statistic':float(ks_2samp(actual,all_depth).statistic),'template_depth_mean':float(actual.mean()),
            'all_controls_depth_mean':float(all_depth.mean())})
        obj.file.close()
    result={'passed':True,'context_specific_effect_shape':list(effects.shape),'outer_target_fold_leakage_check':True,
        'cd4_nested_relative_zero_loss':selection['nested_selection_relative_zero_loss'],
        'v2_cd4_relative_zero_loss':selection['v2_outer_cv_relative_zero_loss'],
        'h1_nested_relative_zero_loss':model['h1_model_selection']['nested_target_holdout_relative_zero_loss'],
        'template_checks':template_checks,'model_sha256':model['model_arrays_sha256'],
        'validator_sha256':build.sha256(Path(__file__)),'scope':'Source-only validation; no challenge perturbation labels or VCC score'}
    build.write_json(A/'model_validation.json',result)
    print('Model axes, source checksums, target folds, weights and all templates validated.',flush=True)
    return manifest,genes,targets,effects


def validate_prediction(manifest,genes,targets,effects):
    ready=json.loads((A/'submission_ready.json').read_text())
    assert ready['official_prep_passed'] and not ready['uploaded']
    assert build.sha256(OUT/'prediction.vcc')==ready['sha256']
    prediction=ad.read_h5ad(OUT/'prediction.h5ad',backed='r')
    assert prediction.shape==(360000,18533)
    np.testing.assert_array_equal(prediction.var_names,genes)
    expected={(c,t):400 for c in ['A','B','C'] for t in targets}
    assert prediction.obs.groupby(['context','target_gene'],observed=True).size().to_dict()==expected
    assert prediction.obs_names.is_unique
    statistics=pd.read_csv(A/'prediction_group_statistics.csv')
    assert len(statistics)==900
    energy=np.sqrt((effects.astype(float)**2).mean(axis=(0,2)))
    coverage=pd.read_csv(A/'target_source_weights.csv').query("context=='A'")
    selected=set(np.argsort(energy)[np.rint(np.linspace(0,299,4)).astype(int)].tolist())
    selected.add(int(np.flatnonzero(coverage.h1_available.to_numpy())[0]))
    selected.add(int(np.flatnonzero(~coverage.cd4_available.to_numpy() & coverage.k562_available.to_numpy())[0]))
    records=[]
    for ci,context in enumerate(['A','B','C']):
        obj=ad.read_h5ad(ROOT/f'context_{context}.h5ad',backed='r')
        template=pd.read_csv(A/f'control_template_{context}.csv')
        sampled=obj.X[template.control_row.to_numpy()].tocsr().astype(float)
        totals=np.asarray(sampled.sum(axis=0)).ravel()
        row_libraries=np.asarray(sampled.sum(axis=1)).ravel()
        stats=pd.read_csv(ROOT/f'analysis/gene_statistics_{context}.csv')
        np.testing.assert_array_equal(stats.gene_name,genes)
        mean_raw=stats.mean_raw_count.to_numpy()
        mean_lib=float(mean_raw.sum())
        for ti in sorted(selected):
            start=(ci*300+ti)*400
            block=prediction.X[start:start+400].tocsr()
            assert np.issubdtype(block.data.dtype,np.integer) and (block.data>0).all()
            factor=np.exp2(effects[ci,ti].astype(float))
            factor[np.flatnonzero(genes==targets[ti])[0]]=.1
            downward=np.minimum(factor,1.)
            rate=(row_libraries.sum()/mean_lib)*mean_raw*np.maximum(factor-1.,0)
            expected_sum=totals*downward+rate
            sum_variance=totals*downward*(1-downward)+rate
            actual=np.asarray(block.sum(axis=0)).ravel()
            # Independent binomial/Poisson expectation, with a conservative integer tail tolerance.
            deviation=np.abs(actual-expected_sum)
            assert (deviation<=8*np.sqrt(sum_variance)+10).all(),(context,targets[ti],float(deviation.max()))
            unchanged=factor==1.
            np.testing.assert_array_equal(actual[unchanged],totals[unchanged])
            entry=statistics[(statistics.context==context)&(statistics.target_gene==targets[ti])].iloc[0]
            assert entry.nnz==block.nnz
            records.append({'context':context,'target_gene':targets[ti],'group_cells':400,'unchanged_gene_sums_exact':True,
                'count_expectations_within_tolerance':True,'gene_tests':len(genes),'mean_library':float(block.sum()/400)})
        obj.file.close()
        print(context,'generated count expectations checked',flush=True)
    prediction.file.close()
    result={'passed':True,'all_groups':900,'cells_per_group':400,'shape':[360000,18533],
        'source_model_validation_passed':True,'official_prep_passed':True,'package_sha256_matches':True,
        'independent_count_expectation_groups':records,'sampled_gene_count_expectation_checks':len(records)*len(genes),
        'validator_sha256':build.sha256(Path(__file__)),'scope':'Numerical validity and simulation checks; no hidden-label accuracy estimate'}
    build.write_json(A/'independent_validation.json',result)
    print('Independent package, annotation and count-expectation checks passed.',flush=True)


def audit_response(genes,targets):
    records=[]
    with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as z:
        effects=z['effects'].copy()
    energy=np.sqrt((effects.astype(float)**2).mean(axis=(0,2)))
    coverage=pd.read_csv(A/'target_source_weights.csv').query("context=='A'")
    selected=set(np.argsort(energy)[np.rint(np.linspace(0,299,4)).astype(int)].tolist())
    selected.add(int(np.flatnonzero(coverage.h1_available.to_numpy())[0]))
    selected.add(int(np.flatnonzero(~coverage.cd4_available.to_numpy() & coverage.k562_available.to_numpy())[0]))
    v2=ad.read_h5ad(ROOT/'experiments/v2_cd4_hybrid/prediction.h5ad',backed='r')
    v3=ad.read_h5ad(OUT/'prediction.h5ad',backed='r')
    for ci,context in enumerate(['A','B','C']):
        control=ad.read_h5ad(ROOT/f'context_{context}.h5ad')
        cpm=normalized_counts(control.X)
        keep=np.flatnonzero(np.asarray(cpm.mean(axis=0)).ravel()>5)
        rank=ControlRankTest(cpm[:,keep])
        for ti in sorted(selected):
            offset=(ci*300+ti)*400
            for name,obj in [('v2',v2),('v3',v3)]:
                block=obj.X[offset:offset+400].tocsr()
                _,q,lfc=rank.test(normalized_counts(block)[:,keep])
                sig=(q<.05)&(genes[keep]!=targets[ti])
                records.append({'context':context,'target_gene':targets[ti],'model':name,'significant_downstream':int(sig.sum()),
                    'up':int((sig&(lfc>0)).sum()),'down':int((sig&(lfc<0)).sum()),'mean_library':float(block.sum()/400)})
            print(context,targets[ti],'v2/v3 downstream calls',[r['significant_downstream'] for r in records[-2:]],flush=True)
        del control,cpm,rank
    v2.file.close(); v3.file.close()
    pd.DataFrame(records).to_csv(A/'v2_v3_response_diagnostics.csv',index=False)
    build.write_json(A/'response_diagnostics_scope.json',{'scope':'Against supplied controls only, NOT response accuracy or a VCC score',
        'targets':targets[sorted(selected)].tolist(),'selection':'Four source-energy quantiles, one H1 target and one K562-only target; deliberately nonrandom',
        'rank_test':'Existing independently checked CPM Mann-Whitney/BH audit; not a local competition scorer'})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model-only',action='store_true')
    parser.add_argument('--audit',action='store_true')
    args=parser.parse_args()
    result=validate_model()
    if not args.model_only:
        validate_prediction(*result)
        if args.audit:
            audit_response(result[1],result[2])
