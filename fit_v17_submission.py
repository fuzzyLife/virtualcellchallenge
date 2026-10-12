"""Freeze V17 production parameters only after supported reserved verification."""
import json
import numpy as np
from threadpoolctl import threadpool_limits
import build_submission as b
from v17_methods import ROOT,A
from run_v17_experiment import load_plan,load_model,write

def main():
 plan=load_plan();selection=json.loads((A/'selection.json').read_text());confirmation=json.loads((A/'confirmation_result.json').read_text());audit=json.loads((A/'confirmation_independent_audit.json').read_text());assert audit['passed'] and confirmation['new_candidate_passes'] and confirmation['candidate']==selection['candidate'] and selection['candidate'] in plan['eligible_finalists'];assert not (A/'reference_model.json').exists()
 model=load_model('all');cv2=model['cv2']*(.5 if selection['candidate']=='learned_half' else 1.);manifest,genes,targets=b.read_spec(ROOT);parent=ROOT/'experiments/v16_control_preserving/analysis';old=json.loads((parent/'reference_model.json').read_text());assert old['model_arrays_sha256']==b.sha256(parent/'reference_model_arrays.npz')
 with np.load(parent/'reference_model_arrays.npz',allow_pickle=False) as z:
  np.testing.assert_array_equal(genes,z['genes']);np.testing.assert_array_equal(targets,z['targets']);effects=z['effects'].copy()
 np.savez_compressed(A/'reference_model_arrays.npz',genes=np.asarray(genes,dtype=str),targets=np.asarray(targets,dtype=str),effects=effects,cv2=cv2)
 deps=set(plan['dependencies'])|{'fit_v17_submission.py','build_v17_submission.py','package_v17_submission.py','validate_v17_submission.py','audit_v17_experiment.py','prepare_v17_confirmation.py','run_v17_confirmation.py'}
 for r in old['context_original_inputs']:assert b.sha256(ROOT/r['file'])==r['sha256']
 write(A/'reference_model.json',{'model':'V17 target-specific response-strength heterogeneity','candidate':selection['candidate'],'response':'frozen_V16','emitter':selection['candidate'],'variation_cv2':cv2.tolist(),'status':'exploratory_public_algorithm_verification','available_public_guards_passed':True,'confirmed_improvement':all(r['complete_primary_metrics'] for r in confirmation['comparisons']),'context_original_inputs':old['context_original_inputs'],'model_arrays_sha256':b.sha256(A/'reference_model_arrays.npz'),'parent_model_arrays_sha256':b.sha256(parent/'reference_model_arrays.npz'),'variation_model_sha256':b.sha256(A/'models/all.npz'),'plan_sha256':b.sha256(A/'plan.json'),'verification_plan_sha256':b.sha256(A/'confirmation_plan.json'),'selection_sha256':b.sha256(A/'selection.json'),'confirmation_sha256':b.sha256(A/'confirmation_result.json'),'independent_experiment_audit_sha256':b.sha256(A/'confirmation_independent_audit.json'),'dependencies':{n:b.sha256(ROOT/n) for n in sorted(deps)},'scope':'Fixed V16 effects and libraries; new shared latent strength only on supported targets. Public source biology reused; no hidden leaderboard guarantee.'});print('V17 production model frozen',selection['candidate'],flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
