"""Freeze V16 identities, methods and release policy before any outcomes."""
from pathlib import Path
import json,platform
import scipy,pandas as pd
import numpy as np
import build_submission as b
from v16_methods import ROOT,A,CANDIDATES

def main():
    A.mkdir(parents=True,exist_ok=True);path=A/'plan.json'
    if path.exists():raise FileExistsError('Preserve V16 prospective plan')
    oldpath=ROOT/'experiments/v15_weighted_transport/analysis/plan.json';old=json.loads(oldpath.read_text());contexts=old['contexts']
    for fi,c in enumerate(contexts):
        c['context']=c['family']+'_v16';c['null_cases']=[]
        for seed in old['seeds']:
            ids=np.random.default_rng(seed+fi+1600).choice(c['left'],600,replace=False).tolist();c['null_cases'].append({'seed':seed,'template':ids[:400],'actual':ids[400:]})
    deps=['v16_methods.py','prepare_v16_experiment.py','run_v16_experiment.py','assess_v16_experiment.py','test_v16_transport.py','v15_ensemble_methods.py','v11_methods.py','v13_quality_methods.py','v7_response_methods.py','v8_methods.py','v6_methods.py','v12_methods.py','v12_cpu_metrics.py','fit_v7_response.py','audit_transfer_validation.py','evaluate_atlas_transfer.py','build_submission.py','test_v14_context.py','v14_context_methods.py']
    guards=old['raw_guards']
    record={'version':16,'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'pandas':pd.__version__},'contexts':contexts,'seeds':old['seeds'],'reserved_targets':old['reserved_targets'],'candidates':CANDIDATES,'dependencies':{n:b.sha256(ROOT/n) for n in deps},'source_cache_sha256':old['source_cache_sha256'],'parent_plan_sha256':b.sha256(oldpath),'raw_guards':guards,
        'design':'Nine frozen arms: V13/calibrated/amplitude-matched uniform response x compact/balanced/residual10 generators. 252 development comparisons, 224 new populations and 28 verified reused baseline populations. Exact copied V15 actual/template/reference IDs; V13 count/profile/DE caches reused only after hashes and original same-input replay verified. No resampling-order assumption.',
        'calibration':'18 fixed bins: significant same-direction family votes 0/1/2+, absolute predicted log2 magnitude <=.05/.05-.2/>.2, fraction agreeing sources <=.75/>.75. Cross-family regression of held effects on predictions; equal family/target/condition weight; nonpanel output genes. Ridge .25*total prediction energy/18 toward gain1; gains clip0..1; protect 2+ votes and >=.2 effects; deployment retains .5 V13 + .5 shrunken V13. All current reserved target labels omitted from gain fitting; known-target source measurements remain legitimate predictor inputs.',
        'balanced':'Per sampled control cell, row-normalize original counts times response factors; thin only expected decreases with binomial removals and redistribute exactly that many molecules by multinomial expected increases. No global mean refit, no resampling of unchanged cells. Zero response and no own target returns original counts exactly; cell-level zero expression cannot be activated. Approximate response means are measured, not forced.',
        'residual10':'Initial .25*expected + .65*Poisson(expected) + .1*row-normalized empirical response counts; retains Poisson coefficient .65 before dual fitting. Same global dual moments as V13; not exact null identity.',
        'null':'Two seeds per family; 400 template and 200 disjoint actual control cells, full independent right controls. Compare generated-vs-right calls, variance and detection; actual200 is descriptive only. New generators require <=1.1 compact false calls, improved mean abs log total-variance ratio, detection MAE <=1.1 compact.',
        'release':'Use unchanged V15 available raw metric and worst-family guard thresholds; require >=3 improved metrics including one DE metric. Equal metric-improvement ranks among passing candidates. Missing NMAE disclosed, never imputed. Each calibrated arm must beat the corresponding exactly per-target L1-amplitude-matched uniform arm on >=2 of fidelity/reach/Jaccard. Uniform arms are eligible candidates under the same guard policy; no new candidate is selected after reserved verification.',
        'verification':'One frozen finalist vs V13 on nine V15-reserved target names at two seeds, plus its amplitude-matched uniform control if the finalist is gene-calibrated; previously inspected biology, not blind independent biology. Same raw guard policy. No retuning or private anchored score estimate. Only a supported exploratory finalist is packaged; no automatic upload.',
        'selection_target_reuse':'Development reuses the 14 V15 targets and source biology; current nine reserved names were not evaluated in V15. This guards algorithm choice, not general biological independence.'}
    b.write_json(path,record);print('V16 prospective plan frozen',b.sha256(path),flush=True)
if __name__=='__main__':main()
