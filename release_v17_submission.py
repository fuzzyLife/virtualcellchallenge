"""Release V17 only after the frozen development and verification gates pass."""
import json,subprocess,sys
import build_submission as b
from v17_methods import ROOT,A

def call(*args):
 print('V17 stage:',*args,flush=True);subprocess.run([sys.executable,*args],cwd=ROOT,check=True)
def main():
 for name in ['models_complete.json','null_evaluation.json','development_evaluation.json']:assert json.loads((A/name).read_text())['complete']
 if not (A/'selection.json').exists():call('assess_v17_experiment.py','development')
 call('audit_v17_experiment.py','development');selection=json.loads((A/'selection.json').read_text())
 if not selection['new_candidate_passes']:
  print('V17 did not pass frozen development; retain V16. Reserved verification remains unopened.',flush=True);return
 if not (A/'confirmation_evaluation.json').exists():call('run_v17_confirmation.py','confirmation','--workers','6')
 if not (A/'confirmation_result.json').exists():call('assess_v17_experiment.py','confirmation')
 call('audit_v17_experiment.py','confirmation');confirmation=json.loads((A/'confirmation_result.json').read_text())
 if not confirmation['new_candidate_passes']:
  print('V17 failed reserved verification; retain V16 without reselection.',flush=True);return
 if not (A/'reference_model.json').exists():call('fit_v17_submission.py')
 if not (A/'preflight.json').exists():call('build_v17_submission.py','--preflight','--workers','6')
 if not (A/'submission_ready.json').exists():call('build_v17_submission.py','--workers','6')
 if not (A/'independent_validation.json').exists():call('validate_v17_submission.py')
 print('V17 independently verified and ready to upload. No upload performed.',flush=True)
if __name__=='__main__':main()
