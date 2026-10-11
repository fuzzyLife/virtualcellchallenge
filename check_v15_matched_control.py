"""Exact mean/blend-matched shuffled-target check for a frozen V15 finalist."""
from pathlib import Path
import json,gc
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
import build_submission as build
from v15_ensemble_methods import ROOT,A,predict
from v15_population_methods import distribution_model,emit
from run_v15_experiment import load_plan,model
from fit_v7_response import read_cache,raw_inputs
from v6_methods import moments
from v12_cpu_metrics import de,direction_metrics

def main():
    plan=load_plan();freeze=json.loads((A/'matched_control_plan.json').read_text());assert freeze['script_sha256']==build.sha256(Path(__file__))
    selection=json.loads((A/'selection.json').read_text());candidate=selection['candidate'];response=json.loads((A/'response_selection.json').read_text())['response']
    path=A/'matched_control_result.json'
    if path.exists():raise FileExistsError('Preserve matched-control result')
    if 'quantile' not in candidate:
        build.write_json(path,{'candidate':candidate,'passed':True,'applicable':False,'scope':'No distribution change to test','selection_sha256':build.sha256(A/'selection.json')});return
    data=read_cache();genes=data['genes'];targets=data['targets'];records=[];shuffle=np.random.default_rng(2026101523).permutation(300)
    amount='quantile10' if candidate.endswith('10') else 'quantile25'
    for ctx in plan['contexts']:
        f=ctx['family'];base,stack=model(data,f);x,meta,measured,_=raw_inputs(f,genes);depth,p,b,_=moments(x[ctx['left']]);query=p*1e6;query[~measured]=np.nan;bank,_=predict(data,base,stack,query);delta,mask,_=distribution_model(data,[f]);ref=x[ctx['right']]
        for cell in [c for c in ctx['cells'] if c['split']=='development']:
            ti=int(np.flatnonzero(targets==cell['target'])[0]);own=int(base['own'][ti]);template=x[cell['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64)
            effect=bank['v13' if candidate.startswith('v13') else response][ti]
            count,_=emit(template,np.arange(400),effect,own,cell['seed']+ti,amount,delta[shuffle[ti]],mask[shuffle[ti]],(libs,p,b))
            truth=de(x[cell['actual']],ref,genes);metrics=direction_metrics(de(count,ref,genes),truth,genes,own)
            records.append({'family':f,'target':cell['target'],'seed':cell['seed'],**metrics});print('V15 matched control',f,cell['target'],cell['seed'],flush=True)
        del x,delta;gc.collect()
    pd.DataFrame(records).to_csv(A/'matched_control_scores.csv',index=False)
    baseline=pd.read_csv(A/'development_scores.csv');baseline=baseline[baseline.candidate==candidate].groupby('family')[['fidelity','reach','jaccard']].mean()
    control=pd.DataFrame(records).groupby('family')[['fidelity','reach','jaccard']].mean();changes=(baseline-control).mean();passed=int((changes>0).sum())>=2
    build.write_json(path,{'candidate':candidate,'passed':passed,'applicable':True,'mean_and_blend_exactly_matched':True,'changes_candidate_minus_shuffled':changes.to_dict(),
        'generated_control_populations':len(records),'selection_sha256':build.sha256(A/'selection.json'),'scores_sha256':build.sha256(A/'matched_control_scores.csv'),'plan_sha256':build.sha256(A/'matched_control_plan.json'),
        'scope':'Predeclared supplementary negative-control check on reused development biology. No reselecting or retuning after failure; not independent confirmation.'})
if __name__=='__main__':
    with threadpool_limits(limits=2):main()
