"""Run independent V64 development contexts, reusing identical V24 counts."""
import subprocess,sys,time,shutil
import v64_response_magnitude as m
from v57_pie import read,write,sha
m.check();m.parent.check();dest=m.A/'development';dest.mkdir(exist_ok=True);waiting=['HEK293T','CD4_Rest_D1','CD4_Stim48_D2','HCT116','H1'];active={};reuse=[]
while waiting or active:
 while waiting and len(active)<3:
  s=waiting.pop(0)
  for src in (m.ROOT/'experiments/v63_identity_guard/analysis/development').glob(s+'__*__v24.*'):
   if src.suffix not in ('.json','.npz'):continue
   target=dest/src.name
   if not target.exists():shutil.copyfile(src,target)
   assert sha(src)==sha(target);reuse.append(dict(file=src.name,sha256=sha(src)))
  log=(m.OUT/(s+'_development.log')).open('w',encoding='utf-8');p=subprocess.Popen([sys.executable,'-X','utf8',str(m.ROOT/'run_v64_counts.py'),'development','--source',s],stdout=log,stderr=subprocess.STDOUT);active[s]=(p,log);print('V64 started',s,flush=True)
 for s,(p,log) in list(active.items()):
  if p.poll() is not None:log.close();assert p.returncode==0,(s,p.returncode);del active[s];print('V64 completed',s,flush=True)
 if waiting or active:time.sleep(5)
write(m.A/'baseline_reuse.json',dict(source='V63/V62 identical original V24 counts and statistics',files=reuse))
subprocess.run([sys.executable,'-X','utf8',str(m.ROOT/'run_v64_counts.py'),'assess-development'],check=True)
