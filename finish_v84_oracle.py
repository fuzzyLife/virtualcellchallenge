"""Run the fixed, conditional oracle diagnosis with two bounded workers."""
import subprocess,sys
import v84_intervention_operator as m
from v57_pie import read,write,sha
assert not read(m.A/'donor_assessment_0.json')['passes']['program25']
assert not (m.A/'oracle').exists(),'Preserve an existing oracle run rather than overwrite its outcomes'
workers=[]
for chunk in range(2):
 log=(m.OUT/f'oracle_chunk{chunk}.log').open('w',encoding='utf-8')
 proc=subprocess.Popen([sys.executable,'-X','utf8','run_v84_oracle.py','run','--chunk',str(chunk),'--chunks','2'],stdout=log,stderr=subprocess.STDOUT)
 workers.append((proc,log))
 print('V84 oracle worker started',chunk,flush=True)
for proc,log in workers:
 rc=proc.wait();log.close();assert rc==0,rc
for script,args in [('run_v84_oracle.py',['assess']),('audit_v84_oracle.py',[])]:
 subprocess.run([sys.executable,'-X','utf8',script,*args],check=True)
write(m.A/'oracle_complete.json',dict(completed=True,non_deployable=True,primary_failed_before_start=True,protocol_sha256=sha(m.A/'oracle_protocol.json'),runner_sha256=sha(__file__)))
