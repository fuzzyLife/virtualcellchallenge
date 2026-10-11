"""Run frozen V16 gates and release only a supported, audited finalist."""
import json,subprocess,sys,time
import build_submission as b
from v16_methods import ROOT,A

def call(*args):
    print('V16 stage:',*args,flush=True);subprocess.run([sys.executable,*args],cwd=ROOT,check=True)

def verify_stage(stage):
    r=json.loads((A/f'{stage}_evaluation.json').read_text());assert r['complete'] and r['plan_sha256']==b.sha256(A/'plan.json')
    assert r['scores_sha256']==b.sha256(A/f'{stage}_scores.csv') and r['panels_sha256']==b.sha256(A/f'{stage}_panels.csv')

def assess(stage):
    file=A/('selection.json' if stage=='development' else 'confirmation_result.json')
    if not file.exists():call('assess_v16_experiment.py',stage)
    r=json.loads(file.read_text());assert r['plan_sha256']==b.sha256(A/'plan.json') and r['scores_sha256']==b.sha256(A/f'{stage}_panels.csv') and r['null_sha256']==b.sha256(A/'null_summary.csv');return r

def main():
    print('Waiting for started model/null/development stages to finish.',flush=True)
    while not all((A/n).exists() for n in ['models_complete.json','null_evaluation.json','development_evaluation.json']):time.sleep(5)
    for stage in ['null','development']:verify_stage(stage)
    models=json.loads((A/'models_complete.json').read_text());assert models['complete'] and models['scores_sha256']==b.sha256(A/'response_scores.csv')
    selection=assess('development');call('audit_v16_experiment.py')
    if selection['candidate']=='v13':
        print('No V16 development candidate passed; retain V13.',flush=True);call('document_v16_submission.py');return
    if not (A/'confirmation_evaluation.json').exists():call('run_v16_experiment.py','confirmation','--workers','6')
    verify_stage('confirmation');confirmation=assess('confirmation');call('audit_v16_experiment.py')
    if not confirmation['new_candidate_passes']:
        print('V16 finalist failed reserved verification; retain V13.',flush=True);call('document_v16_submission.py');return
    if not (A/'reference_model.json').exists():call('fit_v16_submission.py')
    if not (A/'preflight.json').exists():call('build_v16_submission.py','--preflight','--workers','6')
    if not (A/'submission_ready.json').exists():call('build_v16_submission.py','--workers','6')
    if not (A/'independent_validation.json').exists():call('validate_v16_submission.py')
    call('document_v16_submission.py');print('V16 independently validated package ready; no upload performed.',flush=True)
if __name__=='__main__':main()
