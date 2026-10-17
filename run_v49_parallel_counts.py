"""Parallel context execution; identical frozen per-case generator, splits and scores."""
import argparse,json
from concurrent.futures import ProcessPoolExecutor,as_completed
from threadpoolctl import threadpool_limits
import v49_balanced_population as m
import run_v49_count_grid as core
from run_v35_v36 import contexts,clean
from run_v37_ensemble import aggregate
def worker(pair):
 stage,name=pair;core.m=m;core.contexts=lambda:[c for c in contexts() if c['family']==name]
 def write(path,value):
  if path.name in [stage+'_results.json',stage+'_summary.json']:path=path.with_name(name+'_'+path.name)
  m.write(path,clean(value))
 core.write=write
 with threadpool_limits(limits=2):core.counts(stage)
 return name
def main(stage):
 m.check();names=[c['family'] for c in contexts()];record=m.A/'parallel_execution.json'
 receipt={'adapter_sha256':m.b.sha256(m.ROOT/'run_v49_parallel_counts.py'),'frozen_plan_sha256':m.b.sha256(m.A/'plan.json'),'scope':'Four independent context workers, identical existing numerical functions, seeds and splits. Aggregate output paths isolated by context; main process merges only after every worker succeeds. Completed sequential cases reused by hashes. Execution-only optimization after first count outcomes, no changed method or release thresholds.'}
 if record.exists():assert m.read(record)==receipt
 else:m.write(record,receipt)
 with ProcessPoolExecutor(max_workers=4) as pool:
  for f in as_completed([pool.submit(worker,(stage,n)) for n in names]):print('V49 context complete',stage,f.result(),flush=True)
 records=[];panels=[];arms=None
 for name in names:
  r=m.read(m.A/(name+'_'+stage+'_results.json'));records+=r['records'];panels+=r['panels'];arms=r['arms']
 assert len({(r['family'],r['case_index'],r['seed_offset'],r['arm']) for r in records})==len(records)
 m.write(m.A/(stage+'_results.json'),{'records':records,'panels':panels,'arms':arms,'scope':'Identical frozen raw metric calculations, parallel context execution; unscaled public diagnostics.'})
 summary={a:{k:aggregate(panels,a,k) for k in core.METRICS} for a in arms};m.write(m.A/(stage+'_summary.json'),summary);print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);a=p.parse_args();main(a.stage)






