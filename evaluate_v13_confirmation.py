"""Reporting-only adapter: retain missing NMAE for powerless confirmation panels."""
import argparse,importlib,json
from pathlib import Path
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
from v13_methods import ROOT,A
import build_submission as build

def main(workers):
 choice=json.loads((A/'count_selection.json').read_text(encoding='utf-8'));folder=A/choice['validation_folder']
 name={'v13_methods':'evaluate_v13_submission','v13_robust_methods':'evaluate_v13_robust','v13_quality_methods':'evaluate_v13_quality'}[choice['methods_module']]
 module=importlib.import_module(name);original=module.balanced
 def balanced(frame):
  for field in module.FIELDS:
   if field not in frame:frame[field]=np.nan
  return original(frame)
 module.balanced=balanced
 # Original module and script hash remain unchanged, so cached populations and metrics replay exactly.
 module.main('confirmation',workers)
 path=folder/'confirmation_panels.csv';frame=pd.read_csv(path)
 for field in module.FIELDS:
  if field not in frame:frame[field]=np.nan
 frame.to_csv(path,index=False)
 record=json.loads((folder/'confirmation_evaluation.json').read_text(encoding='utf-8'))
 record.update(panels_sha256=build.sha256(path),reporting_adapter_sha256=build.sha256(Path(__file__)))
 build.write_json(folder/'confirmation_evaluation.json',record)
 build.write_json(A/'confirmation_reporting_correction.json',{'change':'Reporting only: preserve ineligible all-None NMAE as NaN columns. No zero substitution, prediction changes, model retuning, split changes or metric changes.',
  'original_evaluator':name+'.py','original_evaluator_sha256':build.sha256(ROOT/(name+'.py')),
  'adapter_sha256':build.sha256(Path(__file__)),'frozen_selection_sha256':build.sha256(A/'count_selection.json')})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);args=p.parse_args()
 with threadpool_limits(limits=2):main(args.workers)
