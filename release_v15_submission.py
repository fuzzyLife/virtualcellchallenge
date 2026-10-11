"""Schedule frozen V15 stages sequentially and preserve failures for inspection."""
from pathlib import Path
import json,subprocess,sys,time
import build_submission as b
ROOT=b.ROOT;A=ROOT/'experiments/v15_weighted_transport/analysis'

def call(*args):
    print('V15 stage:',*args,flush=True)
    subprocess.run([sys.executable,*args],cwd=ROOT,check=True)

def audit():
    call('replay_v15_baseline.py')
    call('audit_v15_experiment.py')

def main():
    while not (A/'response_selection.json').exists():time.sleep(5)
    call('run_v15_experiment.py','null','--workers','6')
    call('run_v15_experiment.py','development','--workers','6')
    call('assess_v15_experiment.py','development')
    selection=json.loads((A/'selection.json').read_text())
    if selection['candidate']=='v13':
        audit();print('No V15 candidate passed development; retain V13.',flush=True);return
    call('check_v15_matched_control.py')
    if not json.loads((A/'matched_control_result.json').read_text())['passed']:
        audit();print('Matched control rejected V15 finalist; retain V13.',flush=True);return
    call('run_v15_experiment.py','confirmation','--workers','6')
    call('assess_v15_experiment.py','confirmation')
    audit()
    if not json.loads((A/'confirmation_result.json').read_text())['new_candidate_passes']:
        print('Reserved algorithm verification rejected V15 finalist; retain V13.',flush=True);return
    call('fit_v15_submission.py')
    call('build_v15_submission.py','--preflight','--workers','6')
    call('build_v15_submission.py','--workers','6')
    call('validate_v15_submission.py')
    print('V15 validated upload package ready; no upload performed.',flush=True)
if __name__=='__main__':main()
