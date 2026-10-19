"""Run sequential V65 selection and locked count tests after all prediction receipts exist."""
import time,subprocess,sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v65_closed_response/analysis'
names=('H1','HCT116','HEK293T','CD4_Rest_D1','CD4_Stim48_D2')
while not all((A/'predictions'/(n+'.json')).exists() for n in names):time.sleep(3)
time.sleep(2)
subprocess.run([sys.executable,'-X','utf8',str(ROOT/'run_v65_experiment.py'),'means'],check=True)
r=json.loads((A/'mean_results.json').read_text())
if r['selected']:
 subprocess.run([sys.executable,'-X','utf8',str(ROOT/'fit_v65_submission.py'),'freeze'],check=True)
 subprocess.run([sys.executable,'-X','utf8',str(ROOT/'run_v65_experiment.py'),'parallel'],check=True)
else:print('No V65 mean finalist; no count or production candidate promoted.',flush=True)
