"""Fit a separately versioned V13 artifact after reporting count confirmation."""
from pathlib import Path
import json,gc
import numpy as np
from threadpoolctl import threadpool_limits
from v13_methods import ROOT,A,banks
from v11_methods import response_for,SETTINGS
from fit_v7_response import read_cache
from v6_methods import moments
import build_submission as build

def main():
 assert not (A/'reference_model.json').exists()
 selection=json.loads((A/'count_selection.json').read_text(encoding='utf-8'));confirmation=json.loads((A/'confirmation_result.json').read_text(encoding='utf-8'))
 assert selection['candidate']==confirmation['candidate'] and confirmation['selection_sha256']==build.sha256(A/'count_selection.json')
 validation=A/selection['validation_folder'];plan=json.loads((validation/'validation_plan.json').read_text(encoding='utf-8'));name=selection['candidate'];setting=SETTINGS['consensus20_noise65']
 from importlib import import_module
 banks=import_module(selection['methods_module']).banks
 data=read_cache();manifest,genes,targets=build.read_spec(ROOT);np.testing.assert_array_equal(data['genes'],genes);np.testing.assert_array_equal(data['targets'],targets)
 bank,training=banks(data,reserved=plan['reserved_targets']);effects=[];inputs=[];diagnostics=[]
 from build_v11_submission import load as load_v11
 _,_,_,_,baseline=load_v11()
 own=np.array([np.flatnonzero(genes==t)[0] for t in targets])
 for ci,context in enumerate(manifest['contexts']):
  x,_,_,_,_,_=build.load_control(ROOT,context,genes,manifest);query=moments(x)[1]*1e6
  effect=np.array([response_for(bank[name],setting,query,int(own[ti]),ti) for ti in range(300)],dtype=np.float32)
  replay=np.array([response_for(bank['v11'],setting,query,int(own[ti]),ti) for ti in range(300)],dtype=np.float32)
  np.testing.assert_array_equal(replay,baseline[ci]);effects.append(effect)
  diagnostics.append({'context':context,'mean_absolute_effect':float(abs(effect).mean()),'changed_from_v11':int((effect!=baseline[ci]).sum()),'mean_absolute_change':float(abs(effect-baseline[ci]).mean())})
  inputs.append({'context':context,'file':f'context_{context}.h5ad','sha256':build.sha256(ROOT/f'context_{context}.h5ad')});del x;gc.collect();print('Fitted V13',context,name,flush=True)
 np.savez_compressed(A/'reference_model_arrays.npz',genes=np.asarray(genes,dtype=str),targets=np.asarray(targets,dtype=str),effects=np.asarray(effects))
 dependencies=['build_v11_submission.py','v13_methods.py','fit_v13_submission.py','build_v13_submission.py','validate_v13_submission.py','package_v13_memory_safe.py',
 'v13_quality_methods.py','evaluate_v13_quality.py','prepare_v13_quality.py','test_v13_quality.py','v13_robust_methods.py','v11_methods.py','v8_methods.py','v7_response_methods.py','fit_v7_response.py','v6_methods.py','build_submission.py','build_v5_submission.py','evaluate_atlas_transfer.py',
 'evaluate_v13_confirmation.py','run_v13_release.py','audit_v13_splits.py','audit_v12_splits.py','select_v13_submission.py','confirm_v13_submission.py','audit_v13_baseline_replay.py','prepare_v13_validation.py','screen_v13_responses.py','screen_v13_robust.py','evaluate_v13_robust.py','prepare_v13_robust.py','test_v13_robust.py','audit_transfer_validation.py','v12_methods.py','v12_cpu_metrics.py','evaluate_v13_submission.py','test_v13_submission.py','assess_v13_submission.py']
 build.write_json(A/'reference_model.json',{'model':'V13 independently validated response model','methods_module':selection['methods_module'],'validation_folder':selection['validation_folder'],'source_checks_file':selection['source_checks_file'],'candidate':name,'setting':setting,
  'confirmed_improvement':confirmation['confirmed_improvement'],'status':'passed_public_raw_guards' if confirmation['confirmed_improvement'] else 'exploratory',
  'model_arrays_sha256':build.sha256(A/'reference_model_arrays.npz'),'context_original_inputs':inputs,
  'source_cache_sha256':build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),
  'validation_plan_file':str((validation/'validation_plan.json').relative_to(A)),'validation_plan_sha256':build.sha256(validation/'validation_plan.json'),'count_selection_sha256':build.sha256(A/'count_selection.json'),
  'confirmation_result_sha256':build.sha256(A/'confirmation_result.json'),'dependencies':{n:build.sha256(ROOT/n) for n in dependencies},
  'contexts':diagnostics,'training':training,'scope':'Public held-family count validation; private ABC treated labels and anchors unavailable. Artifact status must accompany any performance claim.'})
 print('V13 model frozen:',name,confirmation['confirmed_improvement'],flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
