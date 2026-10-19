"""Sequential V66 mean selection and unchanged paired count test."""
import time,subprocess,sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v66_assay_strength/analysis'
names=('H1','HCT116','HEK293T','CD4_Rest_D1','CD4_Stim48_D2')
while not all((A/'predictions'/(n+'.json')).exists() for n in names):time.sleep(3)
time.sleep(2)
subprocess.run([sys.executable,'-X','utf8',str(ROOT/'run_v66_experiment.py'),'means'],check=True)
r=json.loads((A/'mean_results.json').read_text())
if r['selected']:
 subprocess.run([sys.executable,'-X','utf8',str(ROOT/'fit_v66_submission.py'),'freeze'],check=True)
 # Avoid running more than two count-source processes across the two experiments.
while r['selected'] and not (ROOT/'experiments/v65_closed_response/analysis/count_assessment.json').exists():time.sleep(3)
if r['selected']:subprocess.run([sys.executable,'-X','utf8',str(ROOT/'run_v66_experiment.py'),'parallel'],check=True)
else:print('V66 failed the fixed mean screen; no count or production promotion.',flush=True)
