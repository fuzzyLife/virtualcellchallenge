"""Sequence the frozen V89 experimental release without uploading or pushing."""
from pathlib import Path
import sys,subprocess,time
import v89_projection_release as m
from v57_pie import read,write,sha
ROOT=m.ROOT;SCI=sys.executable

def run(script,*args,official=False):
 command=[str(ROOT/'.competition_env/Scripts/python.exe') if official else SCI,'-X','utf8',str(ROOT/script),*map(str,args)];print('V89 stage',script,*args,flush=True);subprocess.run(command,cwd=ROOT,check=True)

def main():
 deadline=time.time()+7200
 while not all((m.A/f'baseline_reuse_620000_{i}.json').exists() for i in range(4)):
  if time.time()>deadline:raise TimeoutError('V89 count workers incomplete')
  time.sleep(15)
 run('mean_v89_release.py')
 for offset in [0,620000]:run('run_v89_counts.py','assess','--offset',offset)
 run('audit_v89_counts.py','--offset',620000)
 old=read(m.old.A/'count_audit_0.json');assert old['passed'];write(m.A/'count_audit_0.json',dict(old,reused_from_V86=True,current_new_populations=0,current_reused_populations=56,original_source_audit_sha256=sha(m.old.A/'count_audit_0.json')))
 run('decide_v89_release.py');decision=read(m.A/'release_decision.json')
 if not decision['experimental_release_allowed']:print('V89 RELEASE BLOCKED BY FROZEN PUBLIC GATES',flush=True);return
 run('fit_v89_submission.py','freeze');run('v89_projection_release.py','fit-production')
 # The immutable V86 inference computes all controls even though only interval25 is deployed.
 impl=m.implementation();impl.fit('production',True)
 run('audit_v89_production.py');run('fit_v89_submission.py','fit');run('build_v89_submission.py','--preflight');run('build_v89_submission.py','--workers',4);run('package_v89_submission.py',official=True);run('validate_v89_submission.py');write(m.A/'experiment_complete.json',dict(completed=True,new_upload_ready=True,uploaded=False,strong_support=decision['strong_support'],script_sha256=sha(__file__)));print('V89 UPLOAD FILE READY AND INDEPENDENTLY VALIDATED',flush=True)
if __name__=='__main__':main()
