"""Prospective V27 calibration and treatment checks; no post-hoc winner picking."""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from v27_methods import ROOT,A,ARMS,FINALIST,load_plan,write
import build_submission as b

def main(stage,family):
 plan=load_plan();tag=stage+'_'+family;dest=A/f'{tag}_result.json';assert not dest.exists()
 evidence=json.loads((A/f'{tag}_evaluation.json').read_text());assert evidence['complete'];assert evidence['scores_sha256']==b.sha256(A/f'{tag}_scores.csv') and evidence['panels_sha256']==b.sha256(A/f'{tag}_panels.csv')
 scores=pd.read_csv(A/f'{tag}_scores.csv');f=pd.read_csv(A/f'{tag}_panels.csv');checks={};changes={};comparisons={}
 if stage=='null':
  rules=plan['null_policy'];rows=scores[scores.candidate==FINALIST];checks['mean_contract']=bool((rows.mean_reference_L1<=.0005).all())
  for _,r in rows.iterrows():
   prefix=str(r.family)+'_'+str(int(r.seed));checks[prefix+'_DE_calls']=bool(r.null_DE_calls<=max(rules['false_call_floor'],r.real_null_DE_calls+rules['real_null_call_allowance']));checks[prefix+'_detection']=bool(r.detection_MAE<=max(rules['detection_MAE_floor'],rules['detection_real_multiplier']*r.real_detection_MAE+rules['detection_real_allowance']));checks[prefix+'_variance']=bool(r.log_variance_MAE<=max(rules['log_variance_MAE_floor'],rules['log_variance_real_multiplier']*r.real_log_variance_MAE+rules['log_variance_real_allowance']))
 else:
  from assess_v24_experiment import compare
  remap=f.copy();remap.candidate=remap.candidate.replace({'v24':'v16',FINALIST:'cdf25','target_shuffled':'shuffled25'});rec=compare(remap,'cdf25',plan['raw_guards']);checks.update(rec['checks']);changes=rec['changes'];checks['mean_contract']=bool((scores[scores.candidate==FINALIST].mean_reference_L1<=.0005).all());base=f[f.candidate=='v24'].set_index(['family','seed']);row=f[f.candidate==FINALIST].set_index(['family','seed']).loc[base.index];checks['shape']=bool((row.shape_wasserstein-base.shape_wasserstein).groupby('family').mean().mean()<0)
  for arm in ARMS:
   other=f[f.candidate==arm].set_index(['family','seed']).loc[base.index];comparisons[arm]={}
   for k in ['pds','fidelity','reach','jaccard','nmae','common_expression_numerator','shape_wasserstein']:
    valid=np.isfinite(other[k])&np.isfinite(row[k]);comparisons[arm][k]=float((row[k]-other[k])[valid].groupby('family').mean().mean()) if valid.any() else None
   if arm in ['anchor_only','calibrated_anchor','calibrated_relative','target_shuffled']:checks['beats_'+arm]=sum(comparisons[arm][k]>0 for k in ['fidelity','reach','jaccard'])>=2
 write(dest,{'passed':all(checks.values()),'checks':checks,'changes':changes,'comparisons_finalist_minus_arm':comparisons,'plan_sha256':b.sha256(A/'plan.json'),'evaluation_sha256':b.sha256(A/f'{tag}_evaluation.json'),'assessor_sha256':b.sha256(Path(__file__))});print(dest.read_text(),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['null','treatment']);p.add_argument('--family',default='H1');a=p.parse_args();main(a.stage,a.family)
