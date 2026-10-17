"""Run independent V52 audits alongside the fitting workflow."""
from pathlib import Path
import subprocess,sys,time
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v52_joint_mean_rank/analysis'
sources=['H1','HCT116','HEK293T','CD4_Rest_D2','CD4_Stim48_D2']
def wait_ready(predicate):
 start=time.monotonic()
 while not predicate():
  if time.monotonic()-start>2400:raise TimeoutError('Expected V52 artifacts not complete; inspect source workflow logs')
  time.sleep(10)
def run(*args):
 print('Audit workflow:',*args,flush=True);subprocess.run([sys.executable,*args],cwd=ROOT,check=True)
wait_ready(lambda:all((A/'labels'/(s+'.json')).exists() for s in sources))
run('audit_v52_joint_response.py','labels')
run('audit_v52_joint_response.py','boundaries')
run('diagnose_v52_label_noise.py')
wait_ready(lambda:(A/'models/all/complete.json').exists() and (A/'evaluation.json').exists())
run('audit_v52_joint_response.py','artifacts')
run('document_v52_joint_response.py')
print('V52 audit workflow completed',flush=True)
