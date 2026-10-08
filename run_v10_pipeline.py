"""Run or resume the reproducible V10 experiment and upload-package preparation."""
from pathlib import Path
import argparse,subprocess,sys,json
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v10_population_transport/analysis'


def run(script,*arguments):
    print('Running',script,*arguments,flush=True)
    subprocess.run([sys.executable,str(ROOT/script),*arguments],cwd=ROOT,check=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--resume',action='store_true');p.add_argument('--workers',type=int,default=4);args=p.parse_args()
    stages=[('validation_plan.json','prepare_v10_validation.py',[]),
            ('count_calibration_summary.csv','evaluate_v10_submission.py',['--stage','calibration','--workers','3']),
            ('decoder_neutral_manifest.json','benchmark_v10_neutral.py',[]),
            ('count_selection.json','select_v10_submission.py',['--stage','calibration']),
            ('count_confirmation_summary.csv','evaluate_v10_submission.py',['--stage','confirmation','--workers','3']),
            ('confirmation_result.json','select_v10_submission.py',['--stage','confirmation']),
            ('reference_model.json','fit_v10_submission.py',[])]
    run('test_v10_submission.py')
    for checkpoint,script,arguments in stages:
        if (A/checkpoint).exists():
            if not args.resume:raise FileExistsError('Existing frozen V10 stage; use --resume to validate and continue')
            print('Preserving completed stage',checkpoint,flush=True)
        else:run(script,*arguments)
    if not (A/'reference_model_arrays.npz').exists():run('restore_v10_model.py')
    if not (A/'corrected_depth_validation.json').exists():run('audit_v10_validation_depths.py')
    # This checks model/dependency hashes and every validation split before generating.
    run('validate_v10_submission.py','--model-only')
    if not (A/'preflight.json').exists():run('build_v10_submission.py','--preflight','--workers',str(args.workers))
    output=ROOT/'experiments/v10_population_transport'
    if not (output/'prediction.vcc').exists():
        if not (output/'prediction.h5ad').exists():run('generate_v10_matrix.py','--workers',str(args.workers))
        run('package_v10_memory_safe.py')
    run('validate_v10_submission.py')
    run('audit_v10_validation_depths.py','--production-only')
    run('analyze_v10_tradeoffs.py')
    run('audit_v10_provenance.py')
    run('document_v10_submission.py')
    print('V10 upload package is ready; this pipeline does not submit it.',flush=True)

if __name__=='__main__':main()
