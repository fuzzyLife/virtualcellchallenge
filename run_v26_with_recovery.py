"""Run the frozen V26 pipeline with only the documented assessor Path adapter."""
from pathlib import Path
import subprocess,sys
import run_v26_pipeline as frozen
from v26_methods import ROOT,load_plan

def main():
 load_plan();original=frozen.call
 def dispatch(script,*args):
  if script=='assess_v26_experiment.py':subprocess.run([sys.executable,str(ROOT/'assess_v26_recovery.py'),*args],cwd=ROOT,check=True)
  else:original(script,*args)
 frozen.call=dispatch;frozen.main()
if __name__=='__main__':main()
