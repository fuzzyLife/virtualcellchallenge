"""Freeze Stack-calibrated production response arrays after reserved approval."""
import json
from pathlib import Path
import numpy as np
import build_submission as b
from v20_production import ROOT,A,verify_production,PROD_DEPS
from v20_methods import target_effects

def main():
 plan=verify_production();assert not (A/'reference_model.json').exists();manifest,genes,targets=b.read_spec(ROOT);genes=np.asarray(genes,dtype=str);targets=np.asarray(targets,dtype=str);parent=ROOT/'experiments/v16_control_preserving/analysis';old=json.loads((parent/'reference_model.json').read_text());assert b.sha256(parent/'reference_model_arrays.npz')==old['model_arrays_sha256']
 with np.load(parent/'reference_model_arrays.npz',allow_pickle=False) as z:np.testing.assert_array_equal(genes,z['genes']);np.testing.assert_array_equal(targets,z['targets']);effects=z['effects'].copy()
 experts=[];changes=[]
 for ci,context in enumerate(manifest['contexts']):
  file=A/'production_experts'/(context+'.npz');side=json.loads(file.with_suffix('.json').read_text());assert side['sha256']==b.sha256(file) and side['production_plan_sha256']==b.sha256(A/'production_plan.json')
  with np.load(file,allow_pickle=False) as z:expert={k:z[k].copy() for k in z.files}
  assert expert['valid'].all();np.testing.assert_array_equal(expert['genes'],genes);np.testing.assert_array_equal(expert['targets'],targets)
  with np.load(A/'production_inputs'/(context+'_query.npz'),allow_pickle=False) as z:query=z['query_mean'].copy()
  # Unused diagnostic aliases avoid allocating duplicate full matrices. The
  # selected stack_mean output depends only on real lfc, shrink and supported.
  expert.update(shuffled_lfc=expert['lfc'],shuffled_shrink=expert['shrink'],raw=expert['lfc'])
  for ti,target in enumerate(targets):
   own=int(np.flatnonzero(genes==target)[0]);base=effects[ci,ti].copy();updated=target_effects(base,expert,ti,query,own)['stack_mean'];assert updated[own]==base[own] and np.max(abs(updated-base))<=.125001;assert np.array_equal(updated[~expert['supported'][ti]],base[~expert['supported'][ti]]);effects[ci,ti]=updated;changes.append({'context':context,'target_gene':str(target),'changed_coordinates':int(np.count_nonzero(updated!=base)),'maximum_absolute_LFC_correction':float(abs(updated-base).max()),'sum_absolute_LFC_correction':float(abs(updated-base).sum())})
  experts.append({'context':context,'path':str(file.relative_to(ROOT)).replace('\\','/'),'sha256':b.sha256(file),'metadata_sha256':b.sha256(file.with_suffix('.json'))})
 np.savez_compressed(A/'reference_model_arrays.npz',genes=genes,targets=targets,effects=effects)
 import pandas as pd
 pd.DataFrame(changes).to_csv(A/'production_response_changes.csv',index=False);confirmation=json.loads((A/'confirmation_result.json').read_text());dependencies=set(json.loads((A/'plan.json').read_text())['dependencies'])|set(PROD_DEPS)|{'prepare_v20_confirmation_inputs.py','audit_v20_confirmation_inputs.py'}
 b.write_json(A/'reference_model.json',{'model':'Stack-calibrated bounded response ensemble V20','candidate':'stack_mean','emitter':'residual10','status':'exploratory_public_algorithm_verification','available_public_guards_passed':True,'confirmed_improvement':all(r['complete_primary_metrics'] for r in confirmation['comparisons']),'context_original_inputs':plan['context_original_inputs'],'model_arrays_sha256':b.sha256(A/'reference_model_arrays.npz'),'parent_model_arrays_sha256':b.sha256(parent/'reference_model_arrays.npz'),'plan_sha256':b.sha256(A/'plan.json'),'production_plan_sha256':b.sha256(A/'production_plan.json'),'selection_sha256':b.sha256(A/'selection.json'),'confirmation_sha256':b.sha256(A/'confirmation_result.json'),'independent_experiment_audit_sha256':b.sha256(A/'confirmation_independent_audit.json'),'source_input_audit_sha256':b.sha256(A/'confirmation_source_input_independent_audit.json'),'experts':experts,'dependencies':{n:b.sha256(ROOT/n) for n in sorted(dependencies)},'pretraining_overlap_audited':False,'scope':'Available development and once-only reserved raw guards passed; reused public biology and missing NMAE remain limitations. No hidden leaderboard ranking guarantee.'});print('Production V20 responses frozen',sum(r['changed_coordinates']>0 for r in changes),'changed groups',flush=True)
if __name__=='__main__':main()
