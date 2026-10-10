"""Complete-panel and worst-family guards for each bounded V12 development stage."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
import build_submission as build
ROOT=build.ROOT;A=ROOT/'experiments/v12_distribution_calibrated/analysis'
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio']
def assess(folder,stage='development'):
 frame=pd.read_csv(folder/f'{stage}_scores.csv');panels=pd.read_csv(folder/f'{stage}_panels.csv');plan=json.loads((folder/'validation_plan.json').read_text())
 names=sorted(frame.candidate.unique());keys=['family','context','target','seed','candidate','reference_mode'];assert not frame.duplicated(keys).any()
 expected={(ctx['family'],ctx['context'],c['target'],c['seed'],n,m) for ctx in plan['contexts'] for c in ctx['cells'] if c['split']==stage for n in names for m in ['sampled400','full_independent']}
 assert set(map(tuple,frame[keys].itertuples(index=False,name=None)))==expected
 results={n:{'passed':True,'modes':{}} for n in names if n!='v11'}
 for mode in ['sampled400','full_independent']:
  f=panels[panels.reference_mode==mode];summary=f.groupby(['candidate','family','context'])[FIELDS].mean().groupby(['candidate','family']).mean();overall=summary.groupby('candidate').mean()
  complete=f.groupby('candidate')[FIELDS].apply(lambda g:bool(np.isfinite(g).all().all())).to_dict()
  for n,r in results.items():
   row=overall.loc[n];base=overall.loc['v11'];checks={'all_six_metrics_complete':complete[n] and complete['v11'],
    'fidelity':bool(row.fidelity>=base.fidelity-.01),'reach':bool(row.reach>=base.reach-.01),
    'nmae':bool(row.nmae<=base.nmae+.01),'pds':bool(row.pds>=base.pds-.02),
    'expression_improves':bool(row.common_expression_ratio<base.common_expression_ratio),
    'family_fidelity':bool(all(summary.loc[(n,fam),'fidelity']>=summary.loc[('v11',fam),'fidelity']-.03 for fam in summary.index.get_level_values('family').unique())),
    'family_reach':bool(all(summary.loc[(n,fam),'reach']>=summary.loc[('v11',fam),'reach']-.03 for fam in summary.index.get_level_values('family').unique()))}
   r['passed'] &= all(checks.values());r['modes'][mode]={'checks':checks,'difference':{k:float(row[k]-base[k]) for k in FIELDS}}
 power=frame[frame.candidate=='v11'].groupby(['family','context','reference_mode']).agg(cases=('target','size'),nmae_eligible=('nmae_eligible','sum'),zero_response=('reference_calls',lambda s:int((s==0).sum())),max_actual_cells=('actual_cells','max')).reset_index()
 power.to_csv(folder/f'{stage}_power.csv',index=False)
 record={'stage':stage,'results':results,'passing_candidates':[n for n,r in results.items() if r['passed']],
  'scores_sha256':build.sha256(folder/f'{stage}_scores.csv'),'plan_sha256':build.sha256(folder/'validation_plan.json'),
  'scope':'Fixed unanchored six-metric guards and complete panels; expression must improve in both reference modes, worst-family fidelity/reach protected. No hosted-score prediction.'}
 build.write_json(folder/f'{stage}_assessment.json',record);print(json.dumps(record,indent=2),flush=True)
 return record
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--folder',default='');p.add_argument('--stage',default='development');a=p.parse_args();assess(A/a.folder,a.stage)
