"""Run disjoint V70 source jobs with bounded memory and explicit failures."""
import subprocess,sys,time
from pathlib import Path
from v70_count_teacher import ROOT,OUT,SOURCES
jobs=list(SOURCES);active={}
while jobs or active:
 while jobs and len(active)<2:
  name=jobs.pop(0);log=(OUT/(name+'_fit.log')).open('w',encoding='utf-8')
  p=subprocess.Popen([sys.executable,'-X','utf8',str(ROOT/'v70_count_teacher.py'),'fit','--source',name],stdout=log,stderr=subprocess.STDOUT)
  active[name]=(p,log);print('Started',name,flush=True)
 for name,(p,log) in list(active.items()):
  if p.poll() is not None:
   log.close();assert p.returncode==0,(name,p.returncode);del active[name];print('Completed',name,flush=True)
 if jobs or active:time.sleep(2)
for name in ['H1','HCT116','HEK293T','CD4_Rest_D1','CD4_Stim48_D2']:
 subprocess.run([sys.executable,'-X','utf8',str(ROOT/'v70_count_teacher.py'),'prepare','--source',name],check=True)
subprocess.run([sys.executable,'-X','utf8',str(ROOT/'v70_count_teacher.py'),'development'],check=True)
