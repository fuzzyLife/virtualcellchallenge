"""Apply prospective raw guards; report incomplete secondary metric panels explicitly."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
import build_submission as build
from v13_methods import A,ROOT
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio']

def assess(stage):
 frame=pd.read_csv(A/f'{stage}_scores.csv');panels=pd.read_csv(A/f'{stage}_panels.csv');plan=json.loads((A/'validation_plan.json').read_text(encoding='utf-8'))
 names=sorted(frame.candidate.unique());keys=['family','context','target','seed','candidate','reference_mode'];assert not frame.duplicated(keys).any()
 expected={(ctx['family'],ctx['context'],c['target'],c['seed'],n,m) for ctx in plan['contexts'] for c in ctx['cells'] if c['split']==stage for n in names for m in ['sampled400','full_independent']}
 assert set(map(tuple,frame[keys].itertuples(index=False,name=None)))==expected
 result={n:{'passed':True,'modes':{}} for n in names if n!='v11'}
 for mode in ['sampled400','full_independent']:
  f=panels[panels.reference_mode==mode];summary=f.groupby(['candidate','family','context'])[FIELDS].mean().groupby(['candidate','family']).mean();overall=summary.groupby('candidate').mean()
  for name,r in result.items():
   row=overall.loc[name];base=overall.loc['v11'];complete=bool(np.isfinite(f[f.candidate.isin([name,'v11'])][FIELDS]).all().all())
   checks={'primary_complete_or_secondary_missing_reported':complete or mode=='sampled400',
    'fidelity':bool(row.fidelity>=base.fidelity-.01),'reach':bool(row.reach>=base.reach-.01),'nmae':bool(row.nmae<=base.nmae+.01),
    'pds':bool(row.pds>=base.pds-.02),'expression_improves':bool(row.common_expression_ratio<base.common_expression_ratio),
    'family_fidelity':bool(all(summary.loc[(name,fam),'fidelity']>=summary.loc[('v11',fam),'fidelity']-.03 for fam in summary.index.get_level_values('family').unique())),
    'family_reach':bool(all(summary.loc[(name,fam),'reach']>=summary.loc[('v11',fam),'reach']-.03 for fam in summary.index.get_level_values('family').unique()))}
   r['passed'] &= all(checks.values());r['modes'][mode]={'checks':checks,'all_six_panels_complete':complete,'difference':{k:float(row[k]-base[k]) if np.isfinite(row[k]-base[k]) else None for k in FIELDS}}
 missing=panels.loc[panels[FIELDS].isna().any(axis=1),['family','context','seed','candidate','reference_mode']].to_dict('records')
 rec={'stage':stage,'passing_candidates':[n for n,r in result.items() if r['passed']],'results':result,'missing_metric_panels':missing,
  'scores_sha256':build.sha256(A/f'{stage}_scores.csv'),'plan_sha256':build.sha256(A/'validation_plan.json'),
  'scope':'Prespecified V13 raw guards. Full-independent references primary and complete; secondary NMAE missingness is disclosed. These are public development data, not hidden competition performance.'}
 build.write_json(A/f'{stage}_assessment.json',rec);print(json.dumps(rec,indent=2),flush=True)
 return rec
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--stage',choices=['development','confirmation'],required=True);assess(p.parse_args().stage)
