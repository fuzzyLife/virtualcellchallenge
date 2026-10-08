"""Run or resume the reproducible V8 experiment and upload-package preparation."""
from pathlib import Path
import argparse,subprocess,sys,json
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v8_corrected_transfer/analysis'


def run(script,*arguments):
    print('Running',script,*arguments,flush=True)
    subprocess.run([sys.executable,str(ROOT/script),*arguments],cwd=ROOT,check=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--resume',action='store_true');p.add_argument('--workers',type=int,default=4);args=p.parse_args()
    stages=[('validation_plan.json','prepare_v8_validation.py',[]),
            ('profile_selection.json','evaluate_v8_submission.py',['--stage','profile']),
            ('count_selection.json','evaluate_v8_submission.py',['--stage','calibration','--workers','2']),
            ('decoder_neutral_manifest.json','benchmark_v8_neutral.py',[]),
            ('safe_count_selection.json','select_v8_safe_generator.py',[]),
            ('confirmation_result.json','evaluate_v8_submission.py',['--stage','confirmation','--workers','2']),
            ('safe_confirmation_result.json','confirm_v8_safe_selection.py',[]),
            ('reference_model.json','fit_v8_submission.py',[])]
    run('test_v8_submission.py')
    for checkpoint,script,arguments in stages:
        if (A/checkpoint).exists():
            if not args.resume:raise FileExistsError('Existing frozen V8 stage; use --resume to validate and continue')
            print('Preserving completed stage',checkpoint,flush=True)
        else:run(script,*arguments)
    # This checks model/dependency hashes and every validation split before generating.
    run('validate_v8_submission.py','--model-only')
    if not (A/'preflight.json').exists():run('build_v8_submission.py','--preflight','--workers',str(args.workers))
    output=ROOT/'experiments/v8_corrected_transfer'
    if not (output/'prediction.vcc').exists():
        if (output/'prediction.h5ad').exists():run('build_v8_submission.py','--package-only')
        else:run('build_v8_submission.py','--workers',str(args.workers))
    run('validate_v8_submission.py')
    run('document_v8_submission.py')
    print('V8 upload package is ready; this pipeline does not submit it.',flush=True)

if __name__=='__main__':main()
