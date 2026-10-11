"""Reconstruct the original V13 baseline on the frozen V15 sampled controls.

Historical V14 output hashes are not a valid comparator because prepare_v14
sampled unsorted control pools before saving sorted pools; V15 samples those
saved pools. No frozen numerical method, identity or outcome is changed here.
"""
from pathlib import Path
import json,gc
from datetime import datetime,timezone
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
import build_submission as b
from run_v15_experiment import load_plan,A,write
from run_v14_experiment import model_cache,count_digest
from v14_context_methods import predict
from v11_methods import emit,SETTINGS
from fit_v7_response import read_cache,raw_inputs
from v6_methods import moments
from v12_cpu_metrics import de,direction_metrics
from audit_transfer_validation import profile_jackknife_chunked
from evaluate_atlas_transfer import load_atlas_model
ROOT=b.ROOT

def main():
    plan=load_plan();old=json.loads((ROOT/'experiments/v14_context_interactions/analysis/plan.json').read_text())
    scores=pd.read_csv(A/'development_scores.csv');scores=scores[scores.candidate=='v13']
    previous=pd.read_csv(ROOT/'experiments/v14_context_interactions/analysis/development_scores.csv');previous=previous[(previous.candidate=='v13')&(previous.reference_mode=='full_independent')]
    pairs=previous.merge(scores,on=['family','target','seed'],suffixes=('_old','_new'));assert len(pairs)==18
    input_checks=[]
    for prior,current in zip(old['contexts'],plan['contexts']):
        assert prior['family']==current['family'] and prior['left']==current['left'] and prior['right']==current['right']
        for c in prior['cells']:
            if c['split']!='development':continue
            v=next(t for t in current['cells'] if t['target']==c['target'] and t['seed']==c['seed'] and t['split']==c['split'])
            assert c['actual']==v['actual']
            input_checks.append({'family':prior['family'],'target':c['target'],'seed':c['seed'],'actual_ids_match':True,'template_ids_match':c['template']==v['template'],'sampled_reference_ids_match':c['reference']==v['reference']})
    assert all(not v['template_ids_match'] for v in input_checks)
    differences=pairs[['family','target','seed','count_sha256_old','count_sha256_new']].copy()
    for f in ['fidelity','reach','jaccard','nmae']:differences[f+'_difference']=pairs[f+'_new']-pairs[f+'_old']
    differences.to_csv(A/'historical_baseline_differences.csv',index=False)
    correction={'created_utc':datetime.now(timezone.utc).isoformat(),'initial_audit_failed':True,'failed_assertion':'All 18 historical V14/V15 baseline count SHA digests must match','historical_digest_matches':int((pairs.count_sha256_old==pairs.count_sha256_new).sum()),'expected_historical_matches':18,
        'reason':'V14 sampled unsorted left/right pools and saved sorted pools; V15 sampled the saved sorted pools. Whole control halves and actual treatment identities match, but sampled template/reference identities do not. Historical count replay is therefore inapplicable.',
        'resolution':'Replay the unchanged original baseline loader, gating and V11 emitter on exact frozen V15 templates/seeds; verify count digest, profiles and DE scores. Keep historical outputs and all V15 methods/inputs/outcomes/selection unchanged.',
        'input_checks':input_checks,'plan_sha256':b.sha256(A/'plan.json'),'initial_auditor_sha256':b.sha256(ROOT/'audit_v15_experiment_initial.py'),'replayer_sha256':b.sha256(Path(__file__)),'outcome_selection_unchanged':True}
    file=A/'audit_correction.json'
    if file.exists():
        saved=json.loads(file.read_text());assert saved['plan_sha256']==correction['plan_sha256'] and saved['replayer_sha256']==correction['replayer_sha256']
    else:write(file,correction)
    data=read_cache();genes=data['genes'];targets=data['targets'];atlas,_=load_atlas_model();records=[]
    for prior,current in zip(old['contexts'],plan['contexts']):
        family=current['family'];assert (ROOT/f'experiments/v14_context_interactions/analysis/models/{family}.npz').exists()
        original=model_cache(data,family);x,meta,measured,_=raw_inputs(family,genes);population=x[current['left']];_,p,bulk,_=moments(population);query=p*1e6;query[~measured]=np.nan
        bank,_=predict(original,query);ref=x[current['right']];old_targets=set(prior['development_targets']);cases=[c for c in current['cells'] if c['split']=='development' and c['target'] in old_targets]
        def one(c):
            ti=int(np.flatnonzero(targets==c['target'])[0]);own=int(original['own'][ti]);template=x[c['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64)
            count,diag=emit(template,np.arange(400),bank['v13'][ti],own,c['seed']+ti,SETTINGS['consensus20_noise65'],prepared=(libs,p,bulk),atlas=atlas)
            record=scores[(scores.family==family)&(scores.target==c['target'])&(scores.seed==c['seed'])].iloc[0]
            assert count_digest(count)==record.count_sha256
            pp,pj=profile_jackknife_chunked(count)
            path=A/'development_cache'/f'{family}_{ti:03d}_{c["seed"]}.npz'
            with np.load(path,allow_pickle=False) as z:
                np.testing.assert_array_equal(pp,z['v13_pred']);np.testing.assert_array_equal(pj,z['v13_jack'])
            metrics=direction_metrics(de(count,ref,genes),de(x[c['actual']],ref,genes),genes,own)
            for f in ['fidelity','reach','jaccard','nmae']:
                actual=metrics[f];wanted=record[f]
                if not np.isfinite(wanted):assert not np.isfinite(actual)
                else:np.testing.assert_allclose(actual,wanted,atol=1e-14,rtol=1e-14)
            out={'family':family,'target':c['target'],'seed':c['seed'],'count_sha256':record.count_sha256,'original_baseline_count_digest_matches':True,'profiles_exact':True,'DE_scores_match':True}
            print('V15 original-baseline same-input replay passed',family,c['target'],c['seed'],flush=True);return out
        with ThreadPoolExecutor(max_workers=6) as pool:records.extend(pool.map(one,cases))
        del x,population,bank,original;gc.collect()
    assert len(records)==18
    write(A/'same_input_baseline_replay.json',{'passed':True,'cases':len(records),'records':records,'plan_sha256':b.sha256(A/'plan.json'),'scores_sha256':b.sha256(A/'development_scores.csv'),'correction_sha256':b.sha256(A/'audit_correction.json'),'script_sha256':b.sha256(Path(__file__)),'scope':'Exact original-baseline response/gating/count emission on V15 inputs; not an exact match of historical V14 populations, which used different templates.'})
    print('V15 all 18 original-baseline same-input replays passed',flush=True)
if __name__=='__main__':
    with threadpool_limits(limits=2):main()
