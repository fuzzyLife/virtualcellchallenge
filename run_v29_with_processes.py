"""Operational V29 process parallelism; frozen algorithms/seeds stay unchanged."""
from pathlib import Path
import argparse,json,sys,subprocess
from joblib import Parallel,delayed,parallel_config
from threadpoolctl import threadpool_limits
from v29_methods import ROOT,A,load_plan,write
import build_submission as b

def worker_init():
 from evaluate_atlas_transfer import load_atlas_model
 from v24_runtime import require_science
 require_science();load_atlas_model()

def run_one(fn,item):
 with threadpool_limits(limits=2):return fn(item)

class Processes:
 def __init__(self,max_workers=2):self.workers=4
 def __enter__(self):return self
 def __exit__(self,*args):return False
 def map(self,fn,items):
  folder=A/'worker_tmp';folder.mkdir(exist_ok=True)
  with parallel_config(backend='loky',inner_max_num_threads=2):
   return Parallel(n_jobs=self.workers,max_nbytes='8M',mmap_mode='r',temp_folder=str(folder),initializer=worker_init)(delayed(run_one)(fn,item) for item in items)

def main():
 plan=load_plan();p=argparse.ArgumentParser();p.add_argument('stage',choices=['pipeline','null','treatment']);p.add_argument('--family',default='H1');args=p.parse_args()
 record={'adapter_sha256':b.sha256(Path(__file__)),'plan_sha256':b.sha256(A/'plan.json'),'frozen_runner_sha256':b.sha256(ROOT/'run_v29_experiment.py'),'workers':4,'BLAS_threads_per_worker':2,'change':'Independent case execution uses joblib loky processes instead of two Python threads. Frozen algorithms, seeds, population identities, numerical runtime, scoring and release rules unchanged; cached cases preserved. Full-axis audit replays use original emitter.'}
 write(A/'process_execution.json',record)
 if args.stage=='pipeline':
  import run_v29_pipeline as pipeline
  original=pipeline.call
  def call(name,*items):
   if name=='run_v29_experiment.py':return subprocess.run([sys.executable,str(Path(__file__).resolve()),*items],cwd=ROOT,check=True)
   return original(name,*items)
  pipeline.call=call
  with threadpool_limits(limits=2):pipeline.main()
 else:
  import run_v29_experiment as runner
  runner.ThreadPoolExecutor=Processes
  with threadpool_limits(limits=2):
   runner.evaluate(args.stage,args.family)
if __name__=='__main__':main()
