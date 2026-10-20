"""Validate the locked V76 candidate; never reselect after confirmation."""
import json,subprocess,sys,time
from pathlib import Path
import v76_response_dilution as m
from v57_pie import read,write,sha

def run(script,*args):
 subprocess.run([sys.executable,'-X','utf8',str(m.ROOT/script),*args],check=True)

def counts(script,label):
 jobs=[c['family'] for c in m.contexts()];active={}
 while jobs or active:
  while jobs and len(active)<2:
   name=jobs.pop(0);log=(m.OUT/f'{name}_{label}.log').open('w',encoding='utf-8');p=subprocess.Popen([sys.executable,'-X','utf8',str(m.ROOT/script),'counts','--source',name],stdout=log,stderr=subprocess.STDOUT);active[name]=(p,log);print('Started V76',label,name,flush=True)
  for name,(p,log) in list(active.items()):
   if p.poll() is not None:
    log.close();assert p.returncode==0,(label,name,p.returncode);del active[name];print('Completed V76',label,name,flush=True)
  if jobs or active:time.sleep(2)
 run(script,'assess')

if __name__=='__main__':
 m.check();assert read(m.A/'selection.json')['selected']=='mixture075';selection_sha=sha(m.A/'selection.json')
 run('v76_response_dilution.py','confirmation')
 log=(m.OUT/'broad.log').open('w',encoding='utf-8');broad=subprocess.Popen([sys.executable,'-X','utf8',str(m.ROOT/'evaluate_v76_broad.py')],stdout=log,stderr=subprocess.STDOUT)
 counts('run_v76_counts.py','confirmation_counts')
 run('audit_v76_counts.py')
 primary=read(m.A/'count_assessment.json');critical=('profile','worst_study','pds','beats_shuffle')
 assert broad.wait()==0;log.close()
 if all(primary['checks'][k] for k in critical) and read(m.A/'broader_assessment.json')['passes']:
  counts('run_v76_repeat_counts.py','repeat_counts');run('audit_v76_counts.py','--repeat')
 else:print('V76 primary critical gates failed; repeat/production not promoted.',flush=True)
 assert sha(m.A/'selection.json')==selection_sha
 write(m.A/'validation_complete.json',dict(selection_unchanged=True,selected='mixture075',primary_critical_gates_passed=all(primary['checks'][k] for k in critical),broad_passed=read(m.A/'broader_assessment.json')['passes'],repeat_completed=(m.A/'count_assessment_repeat.json').exists(),automatic_upload=False))
 print('V76 validation complete',read(m.A/'validation_complete.json'),flush=True)
