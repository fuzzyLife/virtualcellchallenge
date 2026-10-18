"""Resumable V59 public experiments; bounded scheduling, no automatic upload."""
import argparse,subprocess,sys,concurrent.futures
from v59_state_composition import ROOT,A,SOURCES,OUTERS,GROUP,write,sha

def call(*args):
 print('V59',*args,flush=True);subprocess.run([sys.executable,'-X','utf8',str(ROOT/'v59_state_composition.py'),*args],check=True)
def response():
 if not (A/'plan.json').exists():call('freeze')
 subprocess.run([sys.executable,'-X','utf8',str(ROOT/'test_v59_state_composition.py')],check=True)
 for s in SOURCES:call('prepare','--source',s)
 for outer in OUTERS:
  call('states','--outer',outer)
  with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
   jobs=[pool.submit(call,'mechanisms','--outer',outer,'--source',s) for s in SOURCES]
   for job in jobs:job.result()
  for s in SOURCES:
   if GROUP[s]==outer:call('transfer','--source',s)
 call('evaluate');write(A/'response_complete.json',dict(complete=True,plan_sha256=sha(A/'plan.json'),D3_opened=False))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['response','counts']);args=p.parse_args()
 if args.stage=='response':response()
 else:
  with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
   jobs=[pool.submit(subprocess.run,[sys.executable,'-X','utf8',str(ROOT/'evaluate_v59_counts.py'),'run','--source',s],check=True) for s in SOURCES]
   for job in jobs:job.result()
  subprocess.run([sys.executable,'-X','utf8',str(ROOT/'evaluate_v59_counts.py'),'summarize'],check=True)
