"""Terminal rejection from a completed family's irreversible fixed fidelity guard."""
import json
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd,numpy as np
import build_submission as build
ROOT=build.ROOT;A=ROOT/'experiments/v12_distribution_calibrated/analysis/poisson_stage'
def main():
 assert not (A/'development_evaluation.json').exists();plan=json.loads((A/'validation_plan.json').read_text());panels=pd.read_csv(A/'development_panels.csv');scores=pd.read_csv(A/'development_scores.csv');family='H1';summary=panels[panels.family.eq(family)].groupby(['reference_mode','candidate']).fidelity.mean();checks={}
 expected={(c['target'],c['seed'],n,m) for ctx in plan['contexts'] if ctx['family']==family for c in ctx['cells'] if c['split']=='development' for n in plan['emitters'] for m in ['sampled400','full_independent']}
 sub=scores[scores.family.eq(family)];assert set(map(tuple,sub[['target','seed','candidate','reference_mode']].itertuples(index=False,name=None)))==expected
 tolerance=float(plan['guard_tolerances']['family_fidelity'])
 for name in plan['emitters']:
  if name=='v11':continue
  delta={mode:float(summary.loc[(mode,name)]-summary.loc[(mode,'v11')]) for mode in ['sampled400','full_independent']};assert all(v < -tolerance for v in delta.values());checks[name]={'terminal_rejection':True,'family_fidelity_difference':delta,'maximum_allowed_loss':tolerance}
 build.write_json(A/'futility_stop.json',{'created_utc':datetime.now(timezone.utc).isoformat(),'terminal_rejection':True,'results':checks,'completed_family':family,'completed_family_target_seed_cases':len(sub)//len(plan['emitters'])//2,'partial_count_cache_cases':len(list((A/'count_cache').glob('*_development.json'))),'plan_sha256':build.sha256(A/'validation_plan.json'),'scores_sha256':build.sha256(A/'development_scores.csv'),'panels_sha256':build.sha256(A/'development_panels.csv'),'script_sha256':build.sha256(Path(__file__)),'scope':'Rejected early after both seeds of the entire H1 family fail the fixed worst-family guard. Other families cannot repair that failure. Partial later-family caches are retained; no claim of a completed full development run or overall score. No confirmation opened.'})
 print(json.dumps(checks,indent=2),flush=True)
if __name__=='__main__':main()
