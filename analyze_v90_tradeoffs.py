"""Describe V90's already-selected public tradeoffs; never change the deployed candidate."""
import csv
from pathlib import Path
import numpy as np
from v57_pie import read,write,sha
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v90_interval_diagnostic/analysis'

def main():
 rows=[];per={}
 for seed in [0,620000]:
  result=read(ROOT/f'experiments/v89_v86_projection/analysis/count_assessment_{seed}.json');records=result['records'];base={r['target']:r for r in records if r['arm']=='v66'};pred={r['target']:r for r in records if r['arm']=='interval25'}
  for t in sorted(base):
   b=base[t];p=pred[t];row=dict(seed=seed,target=t,profile_ratio=p['profile_MSE']/b['profile_MSE'])
   for k in ['profile_MSE','fidelity','reach','jaccard','nmae']:row[k+'_change']=None if b.get(k) is None or p.get(k) is None else p[k]-b[k]
   rows.append(row);per.setdefault(t,{})[seed]=row
 with (A/'public_target_tradeoffs.csv').open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 summary={}
 for key,smaller in [('profile_ratio',True),('fidelity_change',False),('jaccard_change',False),('nmae_change',True)]:
  cut=1 if key=='profile_ratio' else 0;summary[key]=dict(both_seed_improvements=sum(all(v[s][key] is not None and (v[s][key]<cut if smaller else v[s][key]>cut) for s in [0,620000]) for v in per.values()),both_seed_regressions=sum(all(v[s][key] is not None and (v[s][key]>cut if smaller else v[s][key]<cut) for s in [0,620000]) for v in per.values()))
 worst=min((r for r in rows if r['seed']==620000),key=lambda r:r['jaccard_change']);write(A/'public_target_tradeoffs.json',dict(targets=len(per),summary=summary,worst_second_seed_JAC_case=worst,scope='Descriptive analysis after candidate selection. Same reused donor; no effect on production method, target routing, thresholds or claimed accuracy.',script_sha256=sha(__file__)));print('V90 public target diagnostic',summary,flush=True)
if __name__=='__main__':main()
