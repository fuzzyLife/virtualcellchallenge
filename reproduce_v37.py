"""Outcome-known numerical replay into experiments/replays; original outputs preserved.

Requires frozen V35/V36 model artifacts, prior public input caches, and the pinned
scientific environment. Each stage runs in a fresh process to isolate adapters.
"""
import argparse,datetime,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
STAGES=['primary_freeze','primary_prepare','primary_mean','primary_null','correct_null','primary_development','primary_select','primary_audit','recalibrated_freeze','recalibrated_seal','recalibrated_null','recalibrated_development','recalibrated_select','recalibrated_audit','noise_freeze','noise_seal','noise_null','noise_gate','noise_development','noise_select','noise_audit','final_select']
BOOTSTRAP=r"""
import sys,importlib
from pathlib import Path
from threadpoolctl import threadpool_limits
import run_v37_ensemble as r
folder=Path(sys.argv[1]);stage=sys.argv[2];r.OUT=folder;r.A=folder/'analysis';r.CACHE=r.A/'cache'
def run():
 if stage=='correct_null':
  import run_v37_null_corrected as c;c.main();return
 if stage=='final_select':
  import select_v37_final as f;f.main();return
 kind,action=stage.split('_',1)
 if kind=='primary':
  if action=='audit':
   import audit_v37_ensemble as a;a.main()
  elif action=='select':r.decide('development')
  elif action in ('null','development','confirmation'):r.count(action)
  else:getattr(r,action)()
  return
 if kind=='noise':
  import run_v37_noise_recalibrated as n;x=n.x
  if action=='freeze':n.freeze();return
  if action=='gate':n.gate();return
 else:
  import run_v37_recalibrated as x
 if action in ('seal','audit'):
  import audit_v37_recalibrated as a
  a.A=x.A;a.NAME=x.NAME;a.main(action=='seal');return
 if action=='select':x.install();r.decide('development');return
 if action=='development' and kind=='noise' and not r.read(x.A/'null_gate.json')['passed']:
  print('Noise combination rejected at null gate; treatment replay skipped');return
 if action=='select' and kind=='noise' and not (x.A/'development_results.json').exists():return
 getattr(x,action)()
with threadpool_limits(limits=2):run()
"""

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--directory',type=Path);p.add_argument('--stage',choices=['all',*STAGES],required=True);args=p.parse_args()
 allowed=(ROOT/'experiments/replays').resolve();folder=(args.directory or allowed/('v37_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S'))).resolve()
 assert allowed in folder.parents and folder!=allowed,'Use a fresh directory strictly below experiments/replays.'
 folder.mkdir(parents=True,exist_ok=True)
 stages=STAGES if args.stage=='all' else [args.stage]
 for stage in stages:
  if stage=='noise_select' and not (folder/'analysis/noise_recalibrated/development_results.json').exists():continue
  print('Outcome-known replay',stage,folder,flush=True)
  subprocess.run([sys.executable,'-c',BOOTSTRAP,str(folder),stage],cwd=ROOT,check=True)
 print('Replay completed without modifying original experiment outputs:',folder,flush=True)

if __name__=='__main__':main()
