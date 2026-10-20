"""Bounded completion/audit of already frozen V79-V81 development trials."""
import subprocess,sys,time
from pathlib import Path
names={79:'v79_count_innovations',80:'v80_hurdle_controls',81:'v81_residual_temperature'}
pending=set(names)
while pending:
 for v in sorted(pending):
  a=Path('experiments')/names[v]/'analysis'
  if (a/'development_0_assessment.json').exists():
   print('Independent replay starts',v,flush=True)
   subprocess.run([sys.executable,'-X','utf8','audit_shape_trials.py','--version',str(v),'--stage','development'],check=True)
   pending.remove(v);print('Independent replay completed',v,flush=True)
 if pending:time.sleep(5)
