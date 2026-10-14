"""Once-only V23 selection with paired held-family and matched-control guards."""
import argparse,json
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from v23_methods import A,load_plan,write
import build_submission as b

def main(stage):
 plan=load_plan();dest=A/('selection.json' if stage=='development' else 'confirmation_result.json')
 if dest.exists():raise FileExistsError('Preserve V23 once-only assessment')
 evidence=json.loads((A/f'{stage}_evaluation.json').read_text());assert evidence['complete'] and evidence['panels_sha256']==b.sha256(A/f'{stage}_panels.csv')
 f=pd.read_csv(A/f'{stage}_panels.csv');base=f[f.candidate=='v16'].set_index(['family','seed']);row=f[f.candidate=='mean_preserved'].set_index(['family','seed']).loc[base.index];guards=plan['raw_guards'];changes={};checks={};improved={}
 for metric in ['pds','fidelity','reach','jaccard','nmae','shape_wasserstein']:
  valid=np.isfinite(base[metric])&np.isfinite(row[metric]);delta=(row[metric]-base[metric])[valid];changes[metric]=float(delta.groupby('family').mean().mean()) if valid.any() else None
  improved[metric]=changes[metric] is not None and (changes[metric]<0 if metric in ['nmae','shape_wasserstein'] else changes[metric]>0)
  if metric in ['pds','fidelity','reach','jaccard']:checks[metric]=changes[metric] is not None and changes[metric]>=-guards[metric+'_drop']
 changes['relative_expression_worsening']=float(((row.common_expression_numerator-base.common_expression_numerator)/np.maximum(abs(base.common_expression_numerator),1e-12)).groupby('family').mean().mean());improved['expression']=changes['relative_expression_worsening']<0
 checks.update(nmae=changes['nmae'] is not None and changes['nmae']<=guards['nmae_worsening'],expression=changes['relative_expression_worsening']<=guards['relative_expression_worsening'],three_metrics_improve=sum(improved[k] for k in ['pds','fidelity','reach','jaccard','nmae','expression'])>=3,DE_metric_improves=any(improved[k] for k in ['fidelity','reach','jaccard']),shape_improves=improved['shape_wasserstein'])
 for metric in ['fidelity','reach']:checks['family_'+metric]=bool(((row[metric]-base[metric]).groupby('family').mean()>=-guards['family_'+metric+'_drop']).all())
 for control in ['scalar_preserved','independent_preserved','shuffled_preserved']:
  other=f[f.candidate==control].set_index(['family','seed']).loc[base.index];checks['beats_'+control]=sum(float((row[k]-other[k]).groupby('family').mean().mean())>0 for k in ['fidelity','reach','jaccard'])>=2
 null=pd.read_csv(A/'null_scores.csv');checks['exact_null_identity']=all(group.count_sha256.nunique()==1 for _,group in null.groupby(['family','seed']))
 passed=all(checks.values());write(dest,{'candidate':'mean_preserved' if passed else 'v16','new_candidate_passes':passed,'checks':checks,'changes':changes,'improvements':improved,'complete_primary_metrics':bool(np.isfinite(base.nmae).all() and np.isfinite(row.nmae).all()),'missing_nmae_panels':[{'family':str(v[0]),'seed':int(v[1])} for v in base.index[~np.isfinite(base.nmae)]],'panels_sha256':b.sha256(A/f'{stage}_panels.csv'),'plan_sha256':b.sha256(A/'plan.json'),'scope':'New bounded factor pilot, not complete causal identification or predicted anchored competition score.'});print(dest.read_text(),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);main(p.parse_args().stage)
