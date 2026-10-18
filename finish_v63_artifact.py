"""Finish package preparation and independent verification after generation completes."""
import subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v63_identity_guard/analysis'
science=Path(sys.executable);official=ROOT/'.competition_env/Scripts/python.exe'
while not (A/'prediction_validation.json').exists():time.sleep(3)
for interpreter,script in [(official,'package_v63_submission.py'),(science,'validate_v63_submission.py')]:
 print('V63 release phase',script,flush=True)
 subprocess.run([str(interpreter),'-X','utf8',str(ROOT/script)],check=True)
print('V63 official package and independent payload/provenance audit complete',flush=True)
