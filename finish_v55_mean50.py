"""Finish authorized V55 tests and package; never upload or push."""
from pathlib import Path
import argparse,json,subprocess,sys,time
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v55_mean50_diagnostic/analysis'
def call(script,*args,python=None):
 print('V55 stage:',script,*args,flush=True)
 subprocess.run([str(python or sys.executable),str(ROOT/script),*args],cwd=ROOT,check=True)
def main(wait):
 while not (A/'evaluation.json').exists():
  if not wait:raise FileNotFoundError('Wait for the paired evaluation to complete')
  time.sleep(10)
 # Completed artifact is written only after all context processes return.
 for attempt in range(6):
  try:d=json.loads((A/'evaluation.json').read_text());assert len(d['records'])==176;break
  except (json.JSONDecodeError,AssertionError):
   if attempt==5:raise
   time.sleep(2)
 call('run_v55_mean50.py','assess')
 call('audit_v55_evaluation.py')
 call('run_v55_mean50.py','fit')
 call('audit_v55_mean50.py','production')
 call('run_v55_mean50.py','preflight')
 call('run_v55_mean50.py','generate')
 path=A/'release_decision.json';assert not path.exists();path.write_bytes((A/'selection.json').read_bytes())
 call('package_v55_mean50.py',python=ROOT/'.competition_env/Scripts/python.exe')
 call('audit_v55_mean50.py','delivery')
 call('document_v55_mean50.py')
 print('V55 diagnostic package ready; no upload performed.',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--wait',action='store_true');a=p.parse_args();main(a.wait)
