"""Fit the frozen V10 calibration winner on all public intervention sources."""
import json,gc
import numpy as np
from threadpoolctl import threadpool_limits
from fit_v7_response import read_cache
from v10_methods import ROOT,A,response_bank
from v8_methods import gated
import build_submission as build


def main():
    manifest,genes,targets=build.read_spec(ROOT);data=read_cache()
    np.testing.assert_array_equal(data['genes'],genes);np.testing.assert_array_equal(data['targets'],targets)
    selection=json.loads((A/'count_selection.json').read_text());confirmation=json.loads((A/'confirmation_result.json').read_text())
    assert confirmation['frozen_calibration_selection_sha256']==build.sha256(A/'count_selection.json')
    if (A/'reference_model.json').exists():raise FileExistsError('V10 production model already frozen')
    effects=[];inputs=[];diagnostics=[]
    shared_bank,shared_record=response_bank(data)
    for ci,context in enumerate(manifest['contexts']):
        x,guides,depth,_,_,statistics=build.load_control(ROOT,context,genes,manifest)
        query=np.asarray(x.multiply((1e6/depth)[:,None]).mean(axis=0)).ravel()
        bank,record=response_bank(data,context=context) if selection['response']=='v4' else (shared_bank,shared_record)
        value=bank[selection['response']]
        own=np.array([np.flatnonzero(genes==t)[0] for t in targets])
        response=np.array([gated(value[ti].copy(),query,int(own[ti])) for ti in range(len(targets))],dtype=np.float32)
        effects.append(response);diagnostics.append({'context':context,'nonzero_effects':int((abs(response)>1e-7).sum()),
            'mean_absolute_effect':float(abs(response).mean()),'training_families':record['training_families']})
        inputs.append({'context':context,'file':f'context_{context}.h5ad','sha256':build.sha256(ROOT/f'context_{context}.h5ad')})
        del x;gc.collect();print('Fitted V10 context',context,selection['response'],flush=True)
    np.savez_compressed(A/'reference_model_arrays.npz',genes=np.asarray(genes,dtype=str),targets=np.asarray(targets,dtype=str),effects=np.asarray(effects))
    dependencies=['v10_methods.py','prepare_v10_validation.py','evaluate_v10_submission.py','select_v10_submission.py',
        'fit_v10_submission.py','build_v10_submission.py','validate_v10_submission.py','benchmark_v10_neutral.py',
        'v8_methods.py','audit_transfer_validation.py','fit_v7_response.py','v7_response_methods.py','v6_methods.py',
        'build_submission.py','build_v5_submission.py','evaluate_atlas_transfer.py']
    build.write_json(A/'reference_model.json',{'model':'V10 pooled measured interventions and population dispersion transport',
        'candidate':selection['candidate'],'response':selection['response'],'generator':selection['generator'],
        'confirmation_improvement':confirmation['confirmed_improvement'],'confirmation_utility_difference':confirmation['utility_difference'],
        'model_arrays_sha256':build.sha256(A/'reference_model_arrays.npz'),'context_original_inputs':inputs,
        'source_cache_sha256':build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),
        'validation_plan_sha256':build.sha256(A/'validation_plan.json'),
        'count_selection_sha256':build.sha256(A/'count_selection.json'),'confirmation_result_sha256':build.sha256(A/'confirmation_result.json'),
        'decoder_neutral_manifest_sha256':build.sha256(A/'decoder_neutral_manifest.json'),
        'v4_source_components_sha256':build.sha256(ROOT/'experiments/v4_lineage_response/analysis/source_components.npz'),
        'dependencies':{f:build.sha256(ROOT/f) for f in dependencies},'contexts':diagnostics,
        'knockdown':.9,'gating':'Control mean CPM/5 for output and target, clipped [0,1]; own effect zero and own knockdown applied only in the count emitter.',
        'source_policy':'Calibration-frozen choice across contexts. Pooled responses use availability-masked measurements, bounded source confidence and empirical Bayes shrinkage after pooling. Dose variant normalizes measured source knockdown to .9, bounded .8-1.8; invariant variant subtracts half the source-wide response. V4 candidates retain original components. No inferred tissue identity.',
        'scope':'Public development selection; no hidden challenge treated labels or leaderboard score. A failed confirmation is reported, not retuned away.'})
    print('V10 model frozen',selection['candidate'],selection['generator'],flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=2):main()
