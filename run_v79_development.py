"""Bounded-memory V79 development jobs; no promotion after a failed criterion."""
import subprocess,sys,time
import v79_count_innovations as m
jobs=[c['family'] for c in m.contexts()];active={}
while jobs or active:
 while jobs and len(active)<1:
  name=jobs.pop(0);log=(m.OUT/(name+'_development.log')).open('w',encoding='utf-8');p=subprocess.Popen([sys.executable,'-X','utf8',str(m.ROOT/'run_v79_variance.py'),'counts','--stage','development','--source',name],stdout=log,stderr=subprocess.STDOUT);active[name]=(p,log);print('Started V79 development',name,flush=True)
 for name,(p,log) in list(active.items()):
  if p.poll() is not None:
   log.close();assert p.returncode==0,(name,p.returncode);del active[name];print('Completed V79 development',name,flush=True)
 if jobs or active:time.sleep(2)
subprocess.run([sys.executable,'-X','utf8',str(m.ROOT/'run_v79_variance.py'),'assess','--stage','development'],check=True)
