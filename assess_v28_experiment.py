"""V28 real-control envelopes and protected fixed-mean treatment checks."""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from v28_methods import ROOT,A,FINALIST,load_plan,write
import build_submission as b

def compare(f,scores,plan,global_stage):
 base=f[f.candidate=='v24'].set_index(['family','seed']);row=f[f.candidate==FINALIST].set_index(['family','seed']).loc[base.index];g=plan['raw_guards'];checks={};changes={};improved={}
 for k in ['pds','fidelity','reach','jaccard','nmae']:
  valid=np.isfinite(base[k])&np.isfinite(row[k]);values=(row[k]-base[k])[valid];v=float(values.groupby('family').mean().mean()) if valid.any() else None;changes[k]=v;improved[k]=v is not None and (v<0 if k=='nmae' else v>0)
  if k!='nmae':checks[k]=v is not None and v>=-g[k+'_drop']
  elif v is not None:checks[k]=v<=g['nmae_worsening']
 changes['relative_expression_worsening']=float(((row.common_expression_numerator-base.common_expression_numerator)/np.maximum(abs(base.common_expression_numerator),1e-12)).groupby('family').mean().mean());improved['expression']=changes['relative_expression_worsening']<0;checks['expression']=changes['relative_expression_worsening']<=g['relative_expression_worsening']
 for k in ['fidelity','reach']:checks['family_'+k]=bool(((row[k]-base[k]).groupby('family').mean()>=-g['family_'+k+'_drop']).all())
 changes['shape_wasserstein']=float((row.shape_wasserstein-base.shape_wasserstein).groupby('family').mean().mean());checks['mean_contract']=bool((scores[scores.candidate==FINALIST].mean_reference_L1<=.0005).all())
 if global_stage:checks['three_metrics_improve']=sum(improved.values())>=3;checks['DE_metric_improves']=any(improved[k] for k in ['fidelity','reach','jaccard']);checks['shape_improves']=changes['shape_wasserstein']<0
 return {'passed':all(checks.values()),'checks':checks,'changes':changes,'improvements':improved,'complete_nmae_panels':bool(np.isfinite(row.nmae).all()),'missing_nmae_panels':[{'family':str(v[0]),'seed':int(v[1])} for v in row.index[~np.isfinite(row.nmae)]]}

def main(stage,key):
 plan=load_plan();tag='development' if stage=='development' else stage+'_'+key;dest=A/f'{tag}_result.json';assert not dest.exists()
 if stage=='development':
  scores=pd.concat([pd.read_csv(A/f'treatment_{k}_scores.csv') for k in ['H1','HCT116','HEK293T','CD4']],ignore_index=True);f=pd.concat([pd.read_csv(A/f'treatment_{k}_panels.csv') for k in ['H1','HCT116','HEK293T','CD4']],ignore_index=True);scores.to_csv(A/'development_scores.csv',index=False);f.to_csv(A/'development_panels.csv',index=False)
 else:
  ev=json.loads((A/f'{tag}_evaluation.json').read_text());assert ev['complete'] and ev['scores_sha256']==b.sha256(A/f'{tag}_scores.csv') and ev['panels_sha256']==b.sha256(A/f'{tag}_panels.csv');scores=pd.read_csv(A/f'{tag}_scores.csv');f=pd.read_csv(A/f'{tag}_panels.csv')
 if stage=='null':
  checks={};rul=plan['null_policy'];rows=scores[scores.candidate==FINALIST];checks['mean_contract']=bool((rows.mean_reference_L1<=.0005).all())
  for _,r in rows.iterrows():
   pr=str(r.family)+'_'+str(int(r.seed));checks[pr+'_DE_calls']=bool(r.null_DE_calls<=max(rul['false_call_floor'],r.real_null_DE_calls+rul['real_null_call_allowance']));checks[pr+'_detection']=bool(r.detection_MAE<=max(rul['detection_MAE_floor'],rul['detection_real_multiplier']*r.real_detection_MAE+rul['detection_real_allowance']));checks[pr+'_variance']=bool(r.log_variance_MAE<=max(rul['log_variance_MAE_floor'],rul['log_variance_real_multiplier']*r.real_log_variance_MAE+rul['log_variance_real_allowance']))
  result={'passed':all(checks.values()),'checks':checks}
 else:result=compare(f,scores,plan,stage=='development')
 result.update({'plan_sha256':b.sha256(A/'plan.json'),'scores_sha256':b.sha256(A/f'{tag}_scores.csv'),'panels_sha256':b.sha256(A/f'{tag}_panels.csv'),'assessor_sha256':b.sha256(Path(__file__))});write(dest,result);print(dest.read_text(),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['null','treatment','development']);p.add_argument('--family',default='H1');args=p.parse_args();main(args.stage,args.family)
