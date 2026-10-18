"""Run the fixed V63 development sources as soon as audited predictions exist."""
from pathlib import Path
import subprocess,sys,time,shutil
import v63_identity_guard as m
from v57_pie import read,write,sha
m.check();m.parent.check();dest=m.A/'development';dest.mkdir(parents=True,exist_ok=True)
write(m.A/'mean_results.json',dict(selected='guarded_assay',selection='Only fixed primary, no outcome-based selection',parent_mean_results_sha256=sha(m.parent.A/'mean_results.json')))
waiting=['H1','HCT116','HEK293T','CD4_Rest_D1','CD4_Stim48_D2'];active={};completed=[];reuse=[]
while waiting or active:
 for name in list(waiting):
  if len(active)>=2:break
  if not (m.A/'predictions'/(name+'.json')).exists():continue
  for src in (m.parent.A/'development').glob(name+'__*__v24.*'):
   if src.suffix not in ('.json','.npz'):continue
   out=dest/src.name
   if not out.exists():shutil.copyfile(src,out)
   assert sha(out)==sha(src);reuse.append(dict(file=src.name,source_sha256=sha(src)))
  log=(m.OUT/(name+'_development.log')).open('w',encoding='utf-8')
  p=subprocess.Popen([sys.executable,'-X','utf8',str(m.ROOT/'run_v63_counts.py'),'development','--source',name],stdout=log,stderr=subprocess.STDOUT)
  active[name]=(p,log);waiting.remove(name);print('Started V63 counts',name,flush=True)
 for name,(p,log) in list(active.items()):
  if p.poll() is not None:
   log.close();assert p.returncode==0,(name,p.returncode);completed.append(name);del active[name];print('Completed V63 counts',name,flush=True)
 if waiting or active:time.sleep(5)
write(m.A/'baseline_reuse.json',dict(source='V62 development baseline counts, same source populations/template seeds/emitter and raw statistics',files=reuse,source_plan_sha256=sha(m.parent.A/'plan.json'),parent_count_evaluator_sha256=sha(m.ROOT/'run_v62_counts.py')))
subprocess.run([sys.executable,'-X','utf8',str(m.ROOT/'run_v63_counts.py'),'assess-development'],check=True)
