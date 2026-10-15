"""Gated V26 execution: finish authorized training/tests and retain failed evidence."""
import json,subprocess,sys
from v26_methods import ROOT,A,FAMILIES,load_plan,write
import build_submission as b

def call(script,*args):subprocess.run([sys.executable,str(ROOT/script),*args],cwd=ROOT,check=True)
def main():
 load_plan()
 for family in FAMILIES:
  call('train_v26_decoder.py','train','--family',family)
  if not (A/f'null_{family}_evaluation.json').exists():call('run_v26_experiment.py','null','--family',family)
  if not (A/f'null_{family}_result.json').exists():call('assess_v26_experiment.py','null','--family',family)
  rec=json.loads((A/f'null_{family}_result.json').read_text())
  if not rec['passed']:
   write(A/'calibration_stop.json',{'stopped':True,'failed_family':family,'reason':'Frozen real-control null calibration failed; do not open treatment outcomes or reserve','failed_checks':[k for k,v in rec['checks'].items() if not v],'calibration_sha256':b.sha256(A/f'null_{family}_result.json'),'plan_sha256':b.sha256(A/'plan.json'),'new_package_generated':False,'reserve_opened':False});print('V26 stopped at null gate; retain V24',family,flush=True);return
 if not (A/'development_evaluation.json').exists():call('run_v26_experiment.py','development')
 if not (A/'development_result.json').exists():call('assess_v26_experiment.py','development')
 print('V26 development finished; inspect single locked decision before reserve',flush=True)
if __name__=='__main__':main()
