"""Fit all-source V11 only after the frozen candidate passes count confirmation."""
import json,gc
import numpy as np
from threadpoolctl import threadpool_limits
from v11_methods import ROOT,A,SETTINGS,response_bank,response_for
from fit_v7_response import read_cache
from v6_methods import moments
import build_submission as build


def main():
 if (A/'reference_model.json').exists():raise FileExistsError('V11 model already frozen')
 selection=json.loads((A/'count_selection.json').read_text());confirmation=json.loads((A/'confirmation_result.json').read_text())
 assert selection['candidate']==confirmation['candidate'] and confirmation['selection_sha256']==build.sha256(A/'count_selection.json')
 if not confirmation['confirmed_improvement']:raise ValueError('Frozen candidate did not confirm; preserve V10 and diagnose instead of claiming improvement')
 assert selection['candidate']!='v10';setting=selection['setting'];assert setting==SETTINGS[selection['candidate']]
 data=read_cache();manifest,genes,targets=build.read_spec(ROOT);np.testing.assert_array_equal(data['genes'],genes);np.testing.assert_array_equal(data['targets'],targets)
 bank,training=response_bank(data);effects=[];inputs=[];diagnostics=[]
 baseline_model=json.loads((ROOT/'experiments/v10_population_transport/analysis/reference_model.json').read_text())
 assert build.sha256(ROOT/'experiments/v10_population_transport/analysis/reference_model_arrays.npz')==baseline_model['model_arrays_sha256']
 with np.load(ROOT/'experiments/v10_population_transport/analysis/reference_model_arrays.npz') as z:baseline=z['effects']
 for ci,context in enumerate(manifest['contexts']):
  x,_,depth,_,_,_=build.load_control(ROOT,context,genes,manifest);query=moments(x)[1]*1e6
  own=np.array([np.flatnonzero(genes==t)[0] for t in targets])
  effect=np.array([response_for(bank,setting,query,int(own[ti]),ti) for ti in range(300)],dtype=np.float32)
  equivalent=np.array([response_for(bank,SETTINGS['v10'],query,int(own[ti]),ti) for ti in range(300)],dtype=np.float32)
  np.testing.assert_array_equal(equivalent,baseline[ci])
  effects.append(effect);diagnostics.append({'context':context,'nonzero_effects':int((effect!=0).sum()),
   'mean_absolute_effect':float(abs(effect).mean()),'changed_effects_from_v10':int((effect!=baseline[ci]).sum()),
   'mean_absolute_change_from_v10':float(abs(effect-baseline[ci]).mean())})
  inputs.append({'context':context,'file':f'context_{context}.h5ad','sha256':build.sha256(ROOT/f'context_{context}.h5ad')})
  del x;gc.collect();print('Fitted V11 context',context,selection['candidate'],flush=True)
 np.savez_compressed(A/'reference_model_arrays.npz',genes=np.asarray(genes,dtype=str),targets=np.asarray(targets,dtype=str),effects=np.asarray(effects))
 dependencies=['v11_methods.py','prepare_v11_validation.py','evaluate_v11_submission.py','select_v11_submission.py',
  'fit_v11_submission.py','build_v11_submission.py','validate_v11_submission.py','benchmark_v11_neutral.py',
  'v8_methods.py','v7_response_methods.py','fit_v7_response.py','v6_methods.py','build_submission.py','build_v5_submission.py',
  'evaluate_atlas_transfer.py','audit_transfer_validation.py','select_v8_safe_generator.py','audit_v11_cell_power.py']
 build.write_json(A/'reference_model.json',{'model':'V11 independent-source consensus coverage and controlled emission',
  'candidate':selection['candidate'],'setting':setting,'confirmation_improvement':True,
  'model_arrays_sha256':build.sha256(A/'reference_model_arrays.npz'),'context_original_inputs':inputs,
  'source_cache_sha256':build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),
  'validation_plan_sha256':build.sha256(A/'validation_plan.json'),'count_selection_sha256':build.sha256(A/'count_selection.json'),
  'confirmation_result_sha256':build.sha256(A/'confirmation_result.json'),'decoder_neutral_manifest_sha256':build.sha256(A/'decoder_neutral_manifest.json'),
  'v4_source_components_sha256':build.sha256(ROOT/'experiments/v4_lineage_response/analysis/source_components.npz'),
  'baseline_v10_model_sha256':build.sha256(ROOT/'experiments/v10_population_transport/analysis/reference_model.json'),
  'dependencies':{name:build.sha256(ROOT/name) for name in dependencies},'contexts':diagnostics,'training':training,
  'knockdown':.9,'scope':'Frozen development choice assessed once on count-confirmation; public source summaries previously inspected. H1 reused biology. No hidden treated labels or guaranteed leaderboard gain.'})
 print('V11 model frozen',selection['candidate'],flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
