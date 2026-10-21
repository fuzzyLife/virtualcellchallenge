"""Bounded, resumable V84 donor test; no candidate changes after scoring."""
import subprocess,sys,time
from pathlib import Path
from v57_pie import read,write,sha
import v84_intervention_operator as m
assert (m.A/'donor_predictions.json').exists()
for offset in [0,620000]:
 workers=[]
 for chunk in range(2):
  path=m.OUT/f'donor_{offset}_chunk{chunk}.log';log=path.open('w',encoding='utf-8');p=subprocess.Popen([sys.executable,'-X','utf8','run_v84_donor.py','counts','--offset',str(offset),'--chunk',str(chunk),'--chunks','2'],stdout=log,stderr=subprocess.STDOUT);workers.append((p,log));print('V84 donor worker',offset,chunk,'started',flush=True)
 for p,log in workers:
  rc=p.wait();log.close();assert rc==0,rc
 subprocess.run([sys.executable,'-X','utf8','run_v84_donor.py','assess','--offset',str(offset)],check=True)
 subprocess.run([sys.executable,'-X','utf8','audit_v84_donor.py','--offset',str(offset)],check=True)
write(m.A/'donor_test_complete.json',dict(technical_seeds=[0,620000],both_audits_passed=True,all_CD4_response_evidence_excluded=True,new_model_uploaded=False,script_sha256=sha(__file__)))
