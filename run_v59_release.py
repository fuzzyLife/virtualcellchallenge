"""Canonical V59 release experiment, including pre-score closure correction."""
import argparse,concurrent.futures,subprocess,sys
from v59_state_composition import *

def call(script,*args):
 print('V59 corrected:',script,*args,flush=True)
 subprocess.run([sys.executable,'-X','utf8',str(ROOT/script),*args],check=True)

def freeze_correction():
 dest=A/'pre_score_correction.json'
 if dest.exists():return
 assert not (A/'states').exists() and not (A/'mechanisms').exists() and not (A/'response_assessment.json').exists()
 write(dest,dict(original_plan_sha256=sha(A/'plan.json'),wrapper_sha256=sha(ROOT/'run_v59_closed_stage.py'),controller_sha256=sha(Path(__file__)),reason='Within-state multiplicative changes can alter total CPM. Close every reconstructed profile to 1e6 before LFC comparison; all arms/controls receive same definition. Original code and protocol preserved.',timing='No state fits, treatment reconstruction outcomes, transfer predictions or count scores existed at correction freeze.',unchanged='Raw populations, genes, state architecture, source exclusions, pseudocounts, shrinkage, shuffle, blend weights, cohorts, seeds, count generator and release criteria.',operational='Read-only in-process memoization of verified raw matrices avoids repeated decompression; no transformation changed.',D3_opened=False))

def response():
 freeze_correction();call('test_v59_state_composition.py');call('run_v59_closed_stage.py','check')
 for s in SOURCES:call('v59_state_composition.py','prepare','--source',s)
 for outer in OUTERS:
  call('run_v59_closed_stage.py','states','--outer',outer)
  with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
   jobs=[pool.submit(call,'run_v59_closed_stage.py','mechanisms','--outer',outer,'--source',s) for s in SOURCES]
   for j in jobs:j.result()
  for s in SOURCES:
   if GROUP[s]==outer:call('run_v59_closed_stage.py','transfer','--source',s)
 call('run_v59_closed_stage.py','evaluate')
 write(A/'response_complete.json',dict(complete=True,plan_sha256=sha(A/'plan.json'),correction_sha256=sha(A/'pre_score_correction.json'),D3_opened=False))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','response','counts']);a=p.parse_args()
 if a.stage=='freeze':freeze_correction()
 elif a.stage=='response':response()
 else:call('run_v59_experiment.py','counts')
