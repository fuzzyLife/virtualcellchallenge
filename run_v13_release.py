"""Finish frozen development, one-time confirmation, generation and verification."""
import json,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v13_response_programs/analysis'
def run(file,*args):
 print('Running',file,*args,flush=True)
 subprocess.run([sys.executable,str(ROOT/file),*args],check=True)
def main():
 assert not (A/'count_selection.json').exists(),'Preserve the frozen selection; resume individual documented steps.'
 while not (A/'quality_stage/development_evaluation.json').exists():time.sleep(20)
 subprocess.run([sys.executable,'-c',"import assess_v13_submission as m; m.A=m.ROOT/'experiments/v13_response_programs/analysis/quality_stage'; m.assess('development')"],check=True,cwd=str(ROOT))
 run('compare_v13_quality_anchors.py')
 run('audit_v13_baseline_replay.py')
 run('select_v13_submission.py')
 selection=json.loads((A/'count_selection.json').read_text(encoding='utf-8'))
 evaluator={'v13_methods':'evaluate_v13_submission.py','v13_robust_methods':'evaluate_v13_robust.py','v13_quality_methods':'evaluate_v13_quality.py'}[selection['methods_module']]
 run('evaluate_v13_confirmation.py','--workers','4')
 run('confirm_v13_submission.py')
 run('fit_v13_submission.py')
 run('build_v13_submission.py','--preflight','--workers','4')
 run('build_v13_submission.py','--workers','4')
 run('validate_v13_submission.py')
 run('summarize_v13_evidence.py')
 run('document_v13_submission.py')
 print('V13 release finished. No upload performed.',flush=True)
if __name__=='__main__':main()
