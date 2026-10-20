"""Run bounded-memory V73 source fitting, then the frozen development screen."""
import subprocess,sys,time
import v73_count_assay as m
jobs=list(m.SOURCES);active={}
while jobs or active:
 while jobs and len(active)<2:
  name=jobs.pop(0);log=(m.OUT/(name+'_fit.log')).open('w',encoding='utf-8');p=subprocess.Popen([sys.executable,'-X','utf8',str(m.ROOT/'v73_count_assay.py'),'fit','--source',name],stdout=log,stderr=subprocess.STDOUT);active[name]=(p,log);print('Started V73',name,flush=True)
 for name,(p,log) in list(active.items()):
  if p.poll() is not None:
   log.close();assert p.returncode==0,(name,p.returncode);del active[name];print('Completed V73',name,flush=True)
 if jobs or active:time.sleep(2)
for stage in ('prepare','development'):
 subprocess.run([sys.executable,'-X','utf8',str(m.ROOT/'v73_count_assay.py'),stage],check=True)
