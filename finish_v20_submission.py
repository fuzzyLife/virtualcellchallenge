"""Finish the verified V20 build; optionally wait for the three expert checkpoints."""
from pathlib import Path
import argparse,json,time
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v20_stack_transfer/analysis'

def complete():
 try:return all(len(json.loads((A/'production_experts'/(c+'.json')).read_text())['records'])==300 for c in ['A','B','C'])
 except (FileNotFoundError,json.JSONDecodeError):return False

def main(workers,wait):
 from threadpoolctl import threadpool_limits
 from fit_v20_submission import main as fit
 from build_v20_submission import generate,package,load
 from validate_v20_submission import main as validate
 from release_v20_submission import main as release
 from write_v20_release_docs import main as document
 deadline=time.monotonic()+3600
 while not complete():
  if not wait:raise RuntimeError('Complete run_v20_submission_inference.py first, or use --wait-inference')
  if time.monotonic()>deadline:raise TimeoutError('Inference checkpoints not complete after one hour; inspect the inference job')
  time.sleep(5)
 inputs=json.loads((A/'production_input_independent_audit.json').read_text());assert inputs['passed'];print('All experts and raw-input audits complete; starting final V20 build',flush=True)
 with threadpool_limits(limits=2):
  if not (A/'reference_model.json').exists():print('Freezing response arrays',flush=True);fit()
  load()
  if not (A/'preflight.json').exists():print('Running 24-group count preflight',flush=True);generate(True,workers)
  if not (A.parent/'prediction.h5ad').exists():print('Generating all 900 integer-count groups',flush=True);generate(False,workers)
  if not (A.parent/'prediction.vcc').exists():print('Official lossless packaging',flush=True);package()
  if not (A/'independent_validation.json').exists():print('Independent full-payload and response audit',flush=True);validate()
  release();document()
 print('V20 build, official packaging, independent audits and handoff complete',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4,choices=range(1,7));p.add_argument('--wait-inference',action='store_true');a=p.parse_args();main(a.workers,a.wait_inference)
