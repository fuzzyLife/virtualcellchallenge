"""Descriptive V55 study/context summary; no release rule or prediction changes."""
from pathlib import Path
import json
import numpy as np
import run_v55_mean50 as r
A=r.A

def main():
 d=r.read(A/'evaluation.json');rows=[]
 for family in sorted({x['family'] for x in d['records']}):
  for arm in r.ARMS:
   row={'family':family,'study':r.key(family),'arm':arm}
   rr=[x for x in d['records'] if x['family']==family and x['arm']==arm]
   for k in r.METRICS:
    if k in ('pds','common_expression_ratio'):
     v=[p[k] for p in d['panels'] if p['family']==family and p['arm']==arm and p[k] is not None]
    else:
     v=[]
     for target in sorted({x['target'] for x in rr}):
      vv=[x[k] for x in rr if x['target']==target and x[k] is not None]
      if vv:v.append(float(np.mean(vv)))
    row[k]=float(np.mean(v)) if v else None
   rows.append(row)
 studies=[]
 for study in sorted({x['study'] for x in rows}):
  item={'study':study}
  for arm in r.ARMS:
   item[arm]={}
   for k in r.METRICS:
    vals=[x[k] for x in rows if x['study']==study and x['arm']==arm and x[k] is not None]
    item[arm][k]=float(np.mean(vals)) if vals else None
  item['delta']={k:item['mean50'][k]-item['v24'][k] if item['mean50'][k] is not None and item['v24'][k] is not None else None for k in r.METRICS}
  studies.append(item)
 r.write(A/'study_summary.json',{'contexts':rows,'studies':studies,'scope':'Post-outcome descriptive summary; same evaluated populations and frozen predictions, no selection or rule changes.','evaluation_sha256':r.sha(A/'evaluation.json'),'script_sha256':r.sha(Path(__file__))})
 for x in studies:print(x['study'],x['delta'],flush=True)
if __name__=='__main__':main()
