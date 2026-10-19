"""Use the unchanged V65 unique-target count evaluator for V66."""
import argparse,subprocess,sys,time
from threadpoolctl import threadpool_limits
import v66_assay_strength as m
import run_v65_experiment as runner
runner.m=m

def parallel():
 jobs=[c['family'] for c in runner.contexts()];active={}
 while jobs or active:
  while jobs and len(active)<2:
   name=jobs.pop(0);log=(m.OUT/(name+'_counts.log')).open('w',encoding='utf-8');p=subprocess.Popen([sys.executable,'-X','utf8',str(m.ROOT/'run_v66_experiment.py'),'counts','--source',name],stdout=log,stderr=subprocess.STDOUT);active[name]=(p,log);print('Started V66',name,flush=True)
  for name,(p,log) in list(active.items()):
   if p.poll() is not None:log.close();assert p.returncode==0,(name,p.returncode);del active[name];print('Completed V66',name,flush=True)
  if jobs or active:time.sleep(3)
 runner.assess()

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['means','counts','parallel','assess']);p.add_argument('--source');a=p.parse_args()
 with threadpool_limits(limits=2):
  if a.stage=='counts':runner.counts(a.source)
  elif a.stage=='parallel':parallel()
  else:getattr(runner,a.stage)()
