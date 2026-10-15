"""Once-only V26 real-null calibration and frozen treatment release guards."""
import argparse,json
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from v26_methods import ROOT,A,NAMES,load_plan,write
import build_submission as b

def main(stage,family=None):
 plan=load_plan();tag='null_'+family if stage=='null' else stage;dest=A/f'{tag}_result.json';assert not dest.exists();e=json.loads((A/f'{tag}_evaluation.json').read_text());assert e['complete'] and e['scores_sha256']==b.sha256(A/f'{tag}_scores.csv') and e['panels_sha256']==b.sha256(A/f'{tag}_panels.csv');scores=pd.read_csv(A/f'{tag}_scores.csv');f=pd.read_csv(A/f'{tag}_panels.csv');checks={};changes={}
 if stage=='null':
  rules=plan['null_policy'];rows=scores[scores.candidate=='conditional'];checks['mean_contract']=bool((rows.mean_reference_L1<=plan['mean_contract']['maximum_mean_CPM_probability_L1']).all())
  for _,r in rows.iterrows():
   prefix=str(r.family);checks[prefix+'_DE_calls']=bool(r.null_DE_calls<=max(rules['false_call_floor'],r.real_null_DE_calls+rules['real_null_call_allowance']));checks[prefix+'_detection']=bool(r.detection_MAE<=max(rules['detection_MAE_floor'],rules['detection_real_multiplier']*r.real_detection_MAE+rules['detection_real_allowance']));checks[prefix+'_variance']=bool(r.log_variance_MAE<=max(rules['log_variance_MAE_floor'],rules['log_variance_real_multiplier']*r.real_log_variance_MAE+rules['log_variance_real_allowance']))
 else:
  from assess_v24_experiment import compare
  remap=f.copy();remap.candidate=remap.candidate.replace({'v24':'v16','conditional':'cdf25','target_shuffled':'shuffled25'});rec=compare(remap,'cdf25',plan['raw_guards']);checks.update(rec['checks']);changes=rec['changes'];checks['mean_contract']=bool((scores[scores.candidate=='conditional'].mean_reference_L1<=plan['mean_contract']['maximum_mean_CPM_probability_L1']).all());base=f[f.candidate=='v24'].set_index(['family','seed']);row=f[f.candidate=='conditional'].set_index(['family','seed']).loc[base.index];checks['shape']=bool((row.shape_wasserstein-base.shape_wasserstein).groupby('family').mean().mean()<0)
  for control in ['empirical','independent','state_shuffled']:
   other=f[f.candidate==control].set_index(['family','seed']).loc[base.index];checks['beats_'+control]=sum(float((row[k]-other[k]).groupby('family').mean().mean())>0 for k in ['fidelity','reach','jaccard'])>=2
 result={'passed':all(checks.values()),'checks':checks,'changes':changes,'plan_sha256':b.sha256(A/'plan.json'),'evaluation_sha256':b.sha256(A/f'{tag}_evaluation.json'),'assessor_sha256':b.sha256(__file__)};write(dest,result);print(dest.read_text(),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['null','development']);p.add_argument('--family');a=p.parse_args();main(a.stage,a.family)
