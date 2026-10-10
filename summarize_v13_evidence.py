"""Paired family robustness and conditional target-bootstrap uncertainty."""
import json
from pathlib import Path
import numpy as np,pandas as pd
from v13_methods import A,ROOT
import build_submission as build
FIELDS=['fidelity','reach','jaccard','nmae']

def finite_json(value):
 if isinstance(value,dict):return {k:finite_json(v) for k,v in value.items()}
 if isinstance(value,list):return [finite_json(v) for v in value]
 if isinstance(value,(float,np.floating)) and not np.isfinite(value):return None
 return value

def main():
 choice=json.loads((A/'count_selection.json').read_text(encoding='utf-8'));folder=A/choice['validation_folder'];records={};family_rows=[]
 for stage in ['development','confirmation']:
  frame=pd.read_csv(folder/f'{stage}_scores.csv');frame=frame[frame.reference_mode.eq('full_independent')]
  keys=['family','context','target','seed'];base=frame[frame.candidate.eq('v11')].set_index(keys);new=frame[frame.candidate.eq(choice['candidate'])].set_index(keys)
  assert new.index.equals(base.index)
  delta=(new[FIELDS]-base[FIELDS]).reset_index()
  rng=np.random.default_rng(2026101049);samples=[]
  for _ in range(2000):
   panel=[]
   for family,f in delta.groupby('family'):
    names=f.target.unique();weights=dict(zip(names,rng.exponential(size=len(names))))
    for (context,seed),g in f.groupby(['context','seed']):
     value=g[FIELDS].to_numpy(float);w=g.target.map(weights).to_numpy()[:,None];ok=np.isfinite(value)
     weighted=np.divide(np.nansum(value*w,axis=0),(w*ok).sum(axis=0),out=np.full(len(FIELDS),np.nan),where=(w*ok).sum(axis=0)>0)
     panel.append({'family':family,'context':context,'seed':seed,**dict(zip(FIELDS,weighted))})
   d=pd.DataFrame(panel).groupby(['family','context'])[FIELDS].mean().groupby('family').mean().mean()
   samples.append(d.to_numpy())
  samples=np.array(samples);low,high=np.nanquantile(samples,[.025,.975],axis=0)
  panels=pd.read_csv(folder/f'{stage}_panels.csv');panels=panels[panels.reference_mode.eq('full_independent')]
  metrics=FIELDS+['pds','common_expression_ratio'];balanced=panels.groupby(['candidate','family','context'])[metrics].mean().groupby(['candidate','family']).mean()
  diff=balanced.loc[choice['candidate']]-balanced.loc['v11']
  for family,row in diff.iterrows():family_rows.append({'stage':stage,'family':family,**row.to_dict()})
  records[stage]={'conditional_DE_bootstrap_95pct':{field:{'low':float(low[i]),'high':float(high[i])} for i,field in enumerate(FIELDS)},
   'eligible_family_count_per_metric':diff.notna().sum().to_dict(),'point_family_balanced_delta':diff.mean().to_dict(),'family_minimum_delta':diff.min().to_dict(),'family_maximum_delta':diff.max().to_dict(),
   'fixed_families':diff.index.tolist(),'scores_sha256':build.sha256(folder/f'{stage}_scores.csv')}
 pd.DataFrame(family_rows).to_csv(A/'selected_family_robustness.csv',index=False)
 build.write_json(A/'selected_uncertainty.json',{'candidate':choice['candidate'],'results':finite_json(records),'draws':2000,'seed':2026101049,
  'bootstrap':'Bayesian target bootstrap conditional on the available source families; shared target weight across technical strata and seeds. Panel/seed hierarchy follows the evaluation.',
  'limits':'Not a confidence interval over unseen cell lineages. Seeds and cells are not independent biological experiments. PDS/expression are nonlinear panel metrics: only paired family point ranges are reported for them. Selection used development; confirmation was frozen once.',
  'script_sha256':build.sha256(Path(__file__))})
 print('Conditional target uncertainty and family ranges recorded',flush=True)
if __name__=='__main__':main()
