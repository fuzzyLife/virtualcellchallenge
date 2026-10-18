"""Prepare the fixed V61 experimental release, without uploading it."""
import argparse,subprocess,sys,concurrent.futures
import v61_pooled_state_response as v
m=v.m

def call(script,*args,python=None):
 print('V61 release:',script,*args,flush=True);subprocess.run([python or sys.executable,'-X','utf8',str(m.ROOT/script),*args],check=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--experimental',action='store_true');p.add_argument('--workers',type=int,default=4);a=p.parse_args()
 if not (v.A/'production_plan.json').exists():call('fit_v61_submission.py','freeze')
 call('fit_v61_submission.py','check');call('audit_v61_states.py','boundary');call('audit_v61_states.py','final')
 call('v61_pooled_state_response.py','states','--outer','all')
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
  jobs=[pool.submit(call,'v61_pooled_state_response.py','mechanisms','--outer','all','--source',s) for s in m.SOURCES]
  for j in jobs:j.result()
 call('fit_v61_submission.py','fit',*(['--experimental'] if a.experimental else []))
 call('verify_v61_effects.py');call('build_v61_submission.py','--preflight','--workers',str(a.workers));call('build_v61_submission.py','--workers',str(a.workers))
 call('package_v61_submission.py',python=str(m.ROOT/'.competition_env/Scripts/python.exe'));call('validate_v61_submission.py')
 print('V61 upload file and full integrity audit complete. No submission uploaded.',flush=True)
