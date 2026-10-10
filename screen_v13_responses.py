"""Whole-family-held-out summary screen; these are proxies, not hosted scores."""
import json,gc
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
from fit_v7_response import read_cache
from v7_response_methods import FAMILIES,summary_score
from v13_methods import ROOT,A,NAMES,banks
import build_submission as build

def main():
 A.mkdir(parents=True,exist_ok=True)
 data=read_cache();genes=data['genes'];targets=data['targets']
 previous=ROOT/'experiments/v12_distribution_calibrated/analysis/validation_plan.json'
 p=json.loads(previous.read_text(encoding='utf-8'))
 reserved=sorted({int(np.flatnonzero(targets==c['target'])[0]) for ctx in p['contexts'] for c in ctx['cells'] if c['split']=='confirmation'})
 plan={'candidates':NAMES,'reserved_targets':reserved,'previous_plan_sha256':build.sha256(previous),
 'whole_family_exclusion':True,'methods_sha256':build.sha256(ROOT/'v13_methods.py'),'script_sha256':build.sha256(Path(__file__)),
 'selection':'Advance the best two nonbaseline mean family-balanced proxy objectives to actual count validation; summary improvement alone is insufficient.'}
 path=A/'response_screen_plan.json'
 if path.exists():assert json.loads(path.read_text(encoding='utf-8'))==plan
 else:build.write_json(path,plan)
 records=[];details={}
 for family in FAMILIES:
  bank,record=banks(data,str(family),reserved);details[str(family)]=record
  rows=np.flatnonzero((data['row_families']==family)&~np.isin(data['row_targets'],reserved))
  for ri in rows:
   ti=data['row_targets'][ri];mask=data['valid'][ri].copy();own=int(np.flatnonzero(genes==targets[ti])[0]);mask[own]=False
   for name in NAMES:
    metrics=summary_score(bank[name]['consensus20'][ti],data['truth'][ri],data['errors'][ri],data['q'][ri],mask)
    records.append({'family':str(family),'context':int(data['row_contexts'][ri]),'target':str(targets[ti]),'candidate':name,**metrics})
  pd.DataFrame(records).to_csv(A/'response_screen_scores.csv',index=False)
  print('Completed',family,flush=True);del bank;gc.collect()
 df=pd.DataFrame(records);fields=['mse_ratio','fidelity_proxy','jaccard_proxy','sign_reference_DE','selection_objective']
 family=df.groupby(['candidate','family','context'])[fields].mean().groupby(['candidate','family']).mean()
 family.to_csv(A/'response_screen_family.csv');summary=family.groupby('candidate').mean();summary.to_csv(A/'response_screen_summary.csv')
 ranked=summary.drop(index='v11').sort_values('selection_objective',ascending=False)
 build.write_json(A/'response_screen_result.json',{'complete':True,'advance':ranked.index[:2].tolist(),'training':details,
  'plan_sha256':build.sha256(path),'scores_sha256':build.sha256(A/'response_screen_scores.csv'),
  'scope':'Previously viewed public source-summary proxies; no held-family responses in features, no count confirmation outcomes used.'})
 print(summary.to_string(),flush=True)
if __name__=='__main__':
 from pathlib import Path
 with threadpool_limits(limits=2):main()
