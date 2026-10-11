"""Freeze the best supported V15 candidate after reserved algorithm verification."""
from pathlib import Path
import json,gc
import numpy as np
from threadpoolctl import threadpool_limits
import build_submission as build
from v15_ensemble_methods import ROOT,OUT,A,predict
from v15_population_methods import distribution_model
from run_v15_experiment import load_plan,model,write
from fit_v7_response import read_cache
from v6_methods import moments

def main():
    plan=load_plan();selection=json.loads((A/'selection.json').read_text());confirmation=json.loads((A/'confirmation_result.json').read_text())
    assert selection['candidate']!='v13' and confirmation['candidate']==selection['candidate'] and confirmation['new_candidate_passes'],'Retain V13 after failed candidate verification'
    matched=json.loads((A/'matched_control_result.json').read_text());assert matched['passed'] and matched['candidate']==selection['candidate']
    assert not (A/'reference_model.json').exists()
    data=read_cache();manifest,genes,targets=build.read_spec(ROOT);np.testing.assert_array_equal(data['genes'],genes);np.testing.assert_array_equal(data['targets'],targets)
    base,ensemble=model(data,'all');response=json.loads((A/'response_selection.json').read_text())['response'];name=selection['candidate']
    response='v13' if name.startswith('v13') else response;kind='quantile10' if name.endswith('10') else 'quantile25' if name.endswith('25') else 'compact'
    effects=[];inputs=[]
    with np.load(ROOT/'experiments/v13_response_programs/analysis/reference_model_arrays.npz',allow_pickle=False) as z:baseline=z['effects'].copy()
    for ci,context in enumerate(manifest['contexts']):
        x,*_=build.load_control(ROOT,context,genes,manifest);query=moments(x)[1]*1e6;bank,_=predict(data,base,ensemble,query)
        np.testing.assert_array_equal(bank['v13'],baseline[ci]);effects.append(bank[response])
        inputs.append({'context':context,'file':f'context_{context}.h5ad','sha256':build.sha256(ROOT/f'context_{context}.h5ad')});del x;gc.collect();print('V15 production response fitted',context,flush=True)
    np.savez_compressed(A/'reference_model_arrays.npz',genes=np.asarray(genes,dtype=str),targets=np.asarray(targets,dtype=str),effects=np.asarray(effects,dtype=np.float32))
    distribution_sha=None
    if kind!='compact':
        delta,mask,support=distribution_model(data,[]);np.savez_compressed(A/'distribution_model.npz',delta=delta,mask=mask,support=support);distribution_sha=build.sha256(A/'distribution_model.npz')
    files=list(plan['dependencies'])+['fit_v15_submission.py','build_v15_submission.py','package_v15_submission.py','validate_v15_submission.py','check_v15_matched_control.py','audit_v15_experiment.py']
    old=json.loads((ROOT/'experiments/v13_response_programs/analysis/reference_model.json').read_text())
    write(A/'reference_model.json',{'model':'V15 nested weighted response and source-linked marginal transport','candidate':name,'response':response,'emitter':kind,
        'setting':{'response':'consensus20','amplitude':1.,'noise':.65,'bulk_amplitude':1.,'eligible':True},
        'confirmed_improvement':confirmation['complete_six_metric_verification'],'available_public_guards_passed':confirmation['new_candidate_passes'],'status':'exploratory_public_algorithm_verification',
        'model_arrays_sha256':build.sha256(A/'reference_model_arrays.npz'),'distribution_model_sha256':distribution_sha,'context_original_inputs':inputs,
        'independent_experiment_audit_sha256':build.sha256(A/'independent_experiment_audit.json'),'selection_sha256':build.sha256(A/'selection.json'),'matched_control_sha256':build.sha256(A/'matched_control_result.json'),'confirmation_sha256':build.sha256(A/'confirmation_result.json'),'plan_sha256':build.sha256(A/'plan.json'),
        'source_cache_sha256':plan['source_cache_sha256'],'source_distribution_records':{str(p.relative_to(ROOT)):build.sha256(p) for p in (A/'distribution_sources').glob('*.json')},
        'dependencies':{n:build.sha256(ROOT/n) for n in sorted(set(files))},'baseline_v13_record_sha256':build.sha256(ROOT/'experiments/v13_response_programs/analysis/reference_model.json'),
        'ensemble':ensemble,'scope':'Known perturbation transfer to unperturbed ABC controls; public algorithm verification partly reuses biology, with missing NMAE panels. No hidden leaderboard improvement claim.'})
    print('V15 production model frozen:',name,response,kind,flush=True)
if __name__=='__main__':
    with threadpool_limits(limits=2):main()
