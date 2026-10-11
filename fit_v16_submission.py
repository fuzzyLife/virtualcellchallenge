"""Freeze a supported V16 predictor after once-only reserved verification."""
from pathlib import Path
import json,gc
import numpy as np
from threadpoolctl import threadpool_limits
import build_submission as b
from v16_methods import ROOT,OUT,A,CANDIDATES,predict
from run_v16_experiment import load_plan,model,write
from fit_v7_response import read_cache
from v6_methods import moments

def main():
    plan=load_plan();selection=json.loads((A/'selection.json').read_text());confirmation=json.loads((A/'confirmation_result.json').read_text());audit=json.loads((A/'independent_experiment_audit.json').read_text())
    assert selection['candidate']!='v13' and confirmation['candidate']==selection['candidate'] and confirmation['new_candidate_passes'];assert audit['passed'] and audit['exact_reused_v15_baseline_count_profile_DE_cases']==28 and audit['original_baseline_same_input_regenerations']==6 and audit['bitwise_balanced_null_identity_cases']==6
    assert not (A/'reference_model.json').exists();data=read_cache();manifest,genes,targets=b.read_spec(ROOT);np.testing.assert_array_equal(genes,data['genes']);np.testing.assert_array_equal(targets,data['targets']);bank,calibration=model(data,'all');name=selection['candidate'];response,kind=CANDIDATES[name];effects=[];inputs=[]
    with np.load(ROOT/'experiments/v13_response_programs/analysis/reference_model_arrays.npz',allow_pickle=False) as z:baseline=z['effects'].copy()
    for ci,context in enumerate(manifest['contexts']):
        x,*_=b.load_control(ROOT,context,genes,manifest);query=moments(x)[1]*1e6;pred,factors=predict(bank,calibration,query);np.testing.assert_array_equal(pred['v13'],baseline[ci]);effects.append(pred[response]);inputs.append({'context':context,'file':f'context_{context}.h5ad','sha256':b.sha256(ROOT/f'context_{context}.h5ad')});del x;gc.collect();print('V16 production response fitted',context,flush=True)
    np.savez_compressed(A/'reference_model_arrays.npz',genes=np.asarray(genes,dtype=str),targets=np.asarray(targets,dtype=str),effects=np.asarray(effects,dtype=np.float32))
    deps=set(plan['dependencies'])|{'fit_v16_submission.py','build_v16_submission.py','package_v16_submission.py','validate_v16_submission.py','audit_v16_experiment.py'}
    write(A/'reference_model.json',{'model':'V16 verified response/count transport: '+name,'candidate':name,'response':response,'emitter':kind,'calibration':calibration,'confirmed_improvement':selection['complete_six_metric_verification'] and confirmation['complete_six_metric_verification'],'available_public_guards_passed':True,'status':'exploratory_public_algorithm_verification','model_arrays_sha256':b.sha256(A/'reference_model_arrays.npz'),'context_original_inputs':inputs,
        'plan_sha256':b.sha256(A/'plan.json'),'selection_sha256':b.sha256(A/'selection.json'),'confirmation_sha256':b.sha256(A/'confirmation_result.json'),'independent_experiment_audit_sha256':b.sha256(A/'independent_experiment_audit.json'),'source_cache_sha256':plan['source_cache_sha256'],'dependencies':{n:b.sha256(ROOT/n) for n in sorted(deps)},'scope':'Known-target transfer to ABC controls. Public biology is partly reused; missing NMAE panels prevent complete six-metric verification. No private leaderboard guarantee.'})
    print('V16 production model frozen:',name,response,kind,flush=True)
if __name__=='__main__':
    with threadpool_limits(limits=2):main()
