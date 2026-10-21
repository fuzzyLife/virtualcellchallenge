"""Run the user-requested diagnostic release; never upload or push."""
from pathlib import Path
import sys,subprocess
import v90_diagnostic_release as m
from v57_pie import write,sha
ROOT=m.ROOT

def run(script,*args,official=False):
 print('V90 stage',script,*args,flush=True);subprocess.run([str(ROOT/'.competition_env/Scripts/python.exe') if official else sys.executable,'-X','utf8',str(ROOT/script),*map(str,args)],cwd=ROOT,check=True)
if __name__=='__main__':
 run('v90_diagnostic_release.py','freeze');run('fit_v90_submission.py','freeze');run('v90_diagnostic_release.py','fit-production');run('v90_diagnostic_release.py','fit-production','--shuffled');run('audit_v90_production.py');run('fit_v90_submission.py','fit');run('build_v90_submission.py','--preflight');run('build_v90_submission.py','--workers',4);run('package_v90_submission.py',official=True);run('validate_v90_submission.py');write(m.A/'experiment_complete.json',dict(completed=True,diagnostic_only=True,new_upload_ready=True,uploaded=False,driver_sha256=sha(__file__)));print('V90 DIAGNOSTIC UPLOAD FILE READY',flush=True)
