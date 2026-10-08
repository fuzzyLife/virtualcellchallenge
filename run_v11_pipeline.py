"""Run/resume frozen V11 experiments, generate and audit a new upload file."""
from pathlib import Path
import argparse,subprocess,sys,json
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v11_consensus_coverage/analysis'


def run(script,*args):
 print('Running',script,*args,flush=True);subprocess.run([sys.executable,str(ROOT/script),*args],cwd=ROOT,check=True)


def main():
 p=argparse.ArgumentParser();p.add_argument('--resume',action='store_true');p.add_argument('--workers',type=int,default=4);args=p.parse_args()
 stages=[('validation_plan.json','prepare_v11_validation.py',[]),('development_evaluation.json','evaluate_v11_submission.py',['--stage','development','--workers',str(args.workers)]),
  ('decoder_neutral_manifest.json','benchmark_v11_neutral.py',[]),('development_power_validation.json','audit_v11_cell_power.py',['--stage','development']),('count_selection.json','select_v11_submission.py',['--stage','development']),
  ('confirmation_evaluation.json','evaluate_v11_submission.py',['--stage','confirmation','--workers',str(args.workers)]),
  ('confirmation_power_validation.json','audit_v11_cell_power.py',['--stage','confirmation']),('confirmation_result.json','select_v11_submission.py',['--stage','confirmation']),('reference_model.json','fit_v11_submission.py',[])]
 run('test_v11_submission.py')
 for checkpoint,script,options in stages:
  if (A/checkpoint).exists():
   if not args.resume:raise FileExistsError('Existing frozen stage; use --resume')
   print('Preserving',checkpoint,flush=True)
  else:run(script,*options)
 if not (A/'reference_model_arrays.npz').exists():run('restore_v11_model.py')
 run('validate_v11_submission.py','--model-only')
 if not (A/'preflight.json').exists():run('build_v11_submission.py','--preflight','--workers',str(args.workers))
 out=ROOT/'experiments/v11_consensus_coverage'
 if not (out/'prediction.vcc').exists():
  if not (out/'prediction.h5ad').exists():run('generate_v11_matrix.py','--workers',str(args.workers))
  run('package_v11_memory_safe.py')
 run('validate_v11_submission.py');run('document_v11_submission.py')
 print('V11 independently verified upload file is ready; this pipeline does not submit it.',flush=True)
if __name__=='__main__':main()
