"""Run or resume the reproducible V9 experiment and upload-package preparation."""
from pathlib import Path
import argparse,subprocess,sys,json
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v9_activation_preserving/analysis'


def run(script,*arguments):
    print('Running',script,*arguments,flush=True)
    subprocess.run([sys.executable,str(ROOT/script),*arguments],cwd=ROOT,check=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--resume',action='store_true');p.add_argument('--workers',type=int,default=4);args=p.parse_args()
    stages=[('validation_plan.json','prepare_v9_validation.py',[]),
            ('count_calibration_summary.csv','evaluate_v9_submission.py',['--stage','calibration','--workers','3']),
            ('decoder_neutral_manifest.json','benchmark_v9_neutral.py',[]),
            ('count_selection.json','select_v9_generator.py',[]),
            ('count_confirmation_summary.csv','evaluate_v9_submission.py',['--stage','confirmation','--workers','3']),
            ('confirmation_result.json','confirm_v9_selection.py',[]),
            ('reference_model.json','fit_v9_submission.py',[])]
    run('test_v9_submission.py')
    for checkpoint,script,arguments in stages:
        if (A/checkpoint).exists():
            if not args.resume:raise FileExistsError('Existing frozen V9 stage; use --resume to validate and continue')
            print('Preserving completed stage',checkpoint,flush=True)
        else:run(script,*arguments)
    # This checks model/dependency hashes and every validation split before generating.
    run('validate_v9_submission.py','--model-only')
    if not (A/'preflight.json').exists():run('build_v9_submission.py','--preflight','--workers',str(args.workers))
    output=ROOT/'experiments/v9_activation_preserving'
    if not (output/'prediction.vcc').exists():
        if (output/'prediction.h5ad').exists():run('build_v9_submission.py','--package-only')
        else:run('build_v9_submission.py','--workers',str(args.workers))
    run('validate_v9_submission.py')
    run('audit_v9_provenance.py')
    run('document_v9_submission.py')
    print('V9 upload package is ready; this pipeline does not submit it.',flush=True)

if __name__=='__main__':main()
