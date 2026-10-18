"""Complete extra V63 public count diagnostics without relabeling failed development."""
import argparse,subprocess,sys,time
from pathlib import Path
import v63_identity_guard as m
import run_v62_counts as runner
from v57_pie import sha,read,write
m.data=m.parent.data;runner.m=m

def main(source=None):
 m.check();m.parent.check();dev=m.A/'development_assessment.json'
 while not dev.exists():time.sleep(5)
 verdict=read(dev);stage='confirmation' if verdict['passed'] else 'additional'
 if source:
  with __import__('threadpoolctl').threadpool_limits(limits=2):runner.counts(stage,source)
  return
 path=m.A/'additional_count_plan.json'
 plan=dict(stage=stage,development_passed=verdict['passed'],development_assessment_sha256=sha(dev),candidate='guarded_assay',driver_sha256=sha(__file__),sources=['CD4_Stim48_D2','HEK293T','CD4_Rest_D1','HCT116','H1'],maximum_workers=3,cases='Remaining24 fixed V34 target populations, original seed and +620000, V24/guarded/shuffled arms.',interpretation='If development failed, these are explicitly post-failure diagnostic comparisons on already reused public biology, not confirmation and not retroactive release support. No candidate/parameter reselection; no D3 opening.')
 if path.exists():assert read(path)==plan
 else:write(path,plan)
 waiting=list(plan['sources']);active={}
 while waiting or active:
  while waiting and len(active)<3:
   s=waiting.pop(0);log=(m.OUT/(s+'_'+stage+'.log')).open('w',encoding='utf-8');p=subprocess.Popen([sys.executable,'-X','utf8',str(Path(__file__).resolve()),'--source',s],stdout=log,stderr=subprocess.STDOUT);active[s]=(p,log);print('V63',stage,'started',s,flush=True)
  for s,(p,log) in list(active.items()):
   if p.poll() is not None:log.close();assert p.returncode==0,(s,p.returncode);del active[s];print('V63',stage,'complete',s,flush=True)
  if waiting or active:time.sleep(5)
 with __import__('threadpoolctl').threadpool_limits(limits=2):runner.assess(stage)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source');main(p.parse_args().source)
