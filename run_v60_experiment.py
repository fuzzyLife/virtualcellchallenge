"""Frozen V60 response/count workflow, preserving V59 as a separate experiment."""
import argparse,subprocess,sys,concurrent.futures
import v60_depth_adjusted_states as v
m=v.m

def call(*args):
 print('V60',*args,flush=True);subprocess.run([sys.executable,'-X','utf8',str(m.ROOT/'v60_depth_adjusted_states.py'),*args],check=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['response','counts']);a=p.parse_args()
 if a.stage=='response':
  if not (v.A/'plan.json').exists():call('freeze')
  subprocess.run([sys.executable,'-X','utf8',str(m.ROOT/'test_v60_states.py')],check=True)
  for outer in m.OUTERS:
   call('states','--outer',outer)
   with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    jobs=[pool.submit(call,'mechanisms','--outer',outer,'--source',s) for s in m.SOURCES]
    for j in jobs:j.result()
   for s in m.SOURCES:
    if m.GROUP[s]==outer:call('transfer','--source',s)
  call('evaluate');m.write(v.A/'response_complete.json',dict(complete=True,plan_sha256=m.sha(v.A/'plan.json'),D3_opened=False))
 else:
  with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
   jobs=[pool.submit(call,'count','--source',s) for s in m.SOURCES]
   for j in jobs:j.result()
  call('count_summary')
