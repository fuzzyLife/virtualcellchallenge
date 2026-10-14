"""Schedule unchanged V23 family/seed panels in isolated CPU processes."""
from pathlib import Path
import argparse,json,shutil
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np,pandas as pd
pd.options.future.infer_string=False
import build_submission as b
from v23_methods import ROOT,A,load_plan,write

def one(task):
 family,seed,stage=task
 from threadpoolctl import threadpool_limits
 import run_v23_experiment as runner
 plan=load_plan();ctx=next(c for c in plan['contexts'] if c['family']==family);ctx=dict(ctx,cells=[v for v in ctx['cells'] if v['seed']==seed]);folder=A/'parallel'/f'{family}_{seed}_{stage}';folder.mkdir(parents=True,exist_ok=True);shutil.copyfile(A/'plan.json',folder/'plan.json');(folder/'models').mkdir(exist_ok=True)
 shutil.copyfile(A/'models'/f'{family}.npz',folder/'models'/f'{family}.npz')
 # Only scheduling/output location changes; saved full plan bytes remain identical.
 runner.A=folder;runner.load_plan=lambda:dict(plan,contexts=[ctx])
 with threadpool_limits(limits=2):runner.evaluate(stage,1)
 return str(folder)

def main():
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development']);p.add_argument('--processes',type=int,default=6);args=p.parse_args();plan=load_plan();tasks=[(c['family'],seed,args.stage) for c in plan['contexts'] for seed in sorted({v['seed'] for v in c['cells']})]
 record={'plan_sha256':b.sha256(A/'plan.json'),'scheduler_sha256':b.sha256(Path(__file__)),'frozen_runner_sha256':b.sha256(ROOT/'run_v23_experiment.py'),'tasks':tasks,'processes':args.processes,'BLAS_threads_per_process':2,'numerical_changes':False,'reason':'Python-heavy integer repair does not scale across threads; isolated full family/seed panels retain identical numerical code, cells, seeds, comparators and release rules. Initial threaded attempt produced no completed cases and was stopped; no scoring result informed this scheduling change.'}
 write(A/'parallel_execution.json',record);folders=[]
 with ProcessPoolExecutor(max_workers=args.processes) as pool:
  futures=[pool.submit(one,t) for t in tasks]
  for future in as_completed(futures):folders.append(Path(future.result()));print('V23 process panel complete',folders[-1].name,flush=True)
 scores=pd.concat([pd.read_csv(f/f'{args.stage}_scores.csv') for f in folders],ignore_index=True);panels=pd.concat([pd.read_csv(f/f'{args.stage}_panels.csv') for f in folders],ignore_index=True);scores=scores.sort_values(['family','seed','target','candidate']);panels=panels.sort_values(['family','seed','candidate']);scores.to_csv(A/f'{args.stage}_scores.csv',index=False);panels.to_csv(A/f'{args.stage}_panels.csv',index=False);cache=A/(args.stage+'_cache');cache.mkdir(exist_ok=True)
 for folder in folders:
  for source in (folder/(args.stage+'_cache')).iterdir():
   dest=cache/source.name
   if dest.exists():assert b.sha256(dest)==b.sha256(source)
   else:shutil.copyfile(source,dest)
 fields=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio','shape_wasserstein'];summary=panels.groupby(['candidate','family'])[fields].mean().groupby('candidate').mean();summary.to_csv(A/f'{args.stage}_summary.csv');write(A/f'{args.stage}_evaluation.json',{'complete':True,'populations':len(scores),'case_count':len(scores)//6,'plan_sha256':b.sha256(A/'plan.json'),'scores_sha256':b.sha256(A/f'{args.stage}_scores.csv'),'panels_sha256':b.sha256(A/f'{args.stage}_panels.csv'),'scheduler_sha256':b.sha256(Path(__file__)),'scope':'Identical frozen numeric runner; process scheduling only. Public development biology reused.'});print(summary.to_string(),flush=True)
if __name__=='__main__':main()
