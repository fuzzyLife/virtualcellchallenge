"""Finish this release once its streaming generator has passed full validation.

A one-shot dependent workflow, not a scheduled automation. No upload occurs.
"""
import json,time,subprocess,sys
from pathlib import Path
from v57_pie import ROOT,read,write,sha
A=ROOT/'experiments/v61_pooled_state_response/analysis';OUT=A.parent

def main():
 deadline=time.monotonic()+4*3600
 while True:
  if time.monotonic()>deadline:raise TimeoutError('Generator has not completed within four hours; inspect its process/output before resuming')
  if (A/'prediction_validation.json').exists() and (A/'run_parameters.json').exists() and (OUT/'prediction.h5ad').exists():
   try:
    validation=read(A/'prediction_validation.json');run=read(A/'run_parameters.json');assert run['model_sha256']==sha(A/'reference_model.json');break
   except json.JSONDecodeError:pass
  time.sleep(5)
 print('V61 generation validation complete; starting official packaging',flush=True)
 if not (A/'submission_ready.json').exists():
  subprocess.run([str(ROOT/'.competition_env/Scripts/python.exe'),'-X','utf8',str(ROOT/'package_v61_submission.py')],check=True)
 else:
  ready=read(A/'submission_ready.json');assert sha(OUT/'prediction.vcc')==ready['sha256']
 subprocess.run([sys.executable,'-X','utf8',str(ROOT/'validate_v61_submission.py')],check=True)
 result=read(A/'independent_validation.json');assert result['passed']
 write(A/'packaging_driver.json',dict(passed=True,driver_sha256=sha(__file__),package_sha256=result['package_sha256'],uploaded=False))
 print('V61 experimental package and complete independent validation finished',flush=True)
if __name__=='__main__':main()
