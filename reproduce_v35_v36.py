"""Replay V35/V36 into a fresh namespace without overwriting frozen outcomes.

Requires the large original public input artifacts and a compatible scientific
Python environment. This is outcome-known reproduction, not fresh validation.
"""
import argparse,datetime,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--directory',type=Path);parser.add_argument('--stage',choices=['all','freeze','train35','summary35','audit35','train36','summary36','audit36','audit_cache','mean','null','treatment','audit_counts','select','hybrid_freeze','hybrid_run','hybrid_audit','hybrid_select','contrast_freeze','contrast_run','contrast_audit','contrast_select'],default='all');args=parser.parse_args()
 folder=args.directory or ROOT/'experiments/replays'/('v35_v36_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S'))
 folder=folder.resolve();allowed=(ROOT/'experiments/replays').resolve();assert allowed in folder.parents and folder!=allowed,'Replay must be inside experiments/replays, below its root'
 folder.mkdir(parents=True,exist_ok=True)
 bootstrap="import runpy,sys;from pathlib import Path;import v35_methods as m;import v36_methods as p;root=Path(sys.argv.pop(1));m.A=root/'v35/analysis';m.P36=root/'v36/analysis';m.OUT=m.A.parent;p.A=m.P36;script=sys.argv.pop(1);sys.argv[0]=script;runpy.run_path(script,run_name='__main__')"
 stages=['freeze','train35','summary35','audit35','train36','summary36','audit36','audit_cache','mean','null','treatment','audit_counts','select','hybrid_freeze','hybrid_run','hybrid_audit','hybrid_select','contrast_freeze','contrast_run','contrast_audit','contrast_select'] if args.stage=='all' else [args.stage]
 for stage in stages:
  if stage.startswith('hybrid_') or stage.startswith('contrast_'):
   component,action=stage.split('_',1);suffix='hybrid' if component=='hybrid' else 'contrasts';script=('audit_v35_v36_'+suffix+'.py') if action=='audit' else ('run_v35_v36_'+suffix+'.py');command=[] if action=='audit' else [action]
   if component=='contrast' and action=='select':script='select_v35_v36_contrasts_corrected.py';command=[]
  elif stage=='train36':
   if not (folder/'v36/analysis/normalization_cache_amendment.json').exists():subprocess.run([sys.executable,'-c',bootstrap,str(folder),'run_v35_v36_cached.py','verify'],cwd=ROOT,check=True)
   script='run_v35_v36_cached.py';command=['train36']
  elif stage=='audit_cache':script='audit_v35_v36_cache.py';command=[]
  elif stage=='audit_counts':script='audit_v35_v36_counts.py';command=[]
  elif stage.startswith('audit'):script='audit_v35_v36.py';command=[stage[-2:]]
  elif stage=='select':script='select_v35_v36_release_corrected.py';command=['select']
  else:script='run_v35_v36.py';command=[stage]
  print('Replay',stage,'in',folder,flush=True);subprocess.run([sys.executable,'-c',bootstrap,str(folder),script,*command],cwd=ROOT,check=True)
  if stage=='freeze':
   subprocess.run([sys.executable,'-c',bootstrap,str(folder),'select_v35_v36_release.py','freeze'],cwd=ROOT,check=True)
 print('Outcome-known replay finished:',folder,flush=True)
if __name__=='__main__':main()