"""Reproduce one completed pilot in a fresh output directory; preserve originals.

The original large competition/public source inputs are still required. This is
an outcome-known numerical reproduction, never a new blind validation run.
"""
import argparse,json,shutil
from datetime import datetime,timezone
from functools import lru_cache
from threadpoolctl import threadpool_limits
import v32_methods as m
import v32_runtime as runtime

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('version',choices=['32','33','34']);args=p.parse_args();version=args.version;old=m.ROOT/'experiments'/ {'32':'v32_paired_context','33':'v33_tabular_response','34':'v34_residual_flow'}[version]/'analysis';dest=m.ROOT/'experiments/replays'/f'v{version}_{datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")}'/'analysis';assert not dest.exists();dest.mkdir(parents=True);shutil.copyfile(old/'plan.json',dest/'plan.json');m.write(dest/'reproduction.json',{'original_plan_sha256':m.sha(old/'plan.json'),'new_plan_sha256':m.sha(dest/'plan.json'),'original_analysis':str(old),'output_analysis':str(dest),'runner_sha256':m.sha(m.ROOT/'reproduce_v32_v34.py'),'scope':'Same frozen settings/source data in a new namespace. Outcomes already known; numerical reproduction, not new validation.'});m.old.embeddings=lru_cache(maxsize=1)(m.old.embeddings);runtime.install();print('New reproduction output:',dest,flush=True)
 if version=='32':
  r=json.loads((old/'corrected_results.json').read_text());assert m.sha(m.ROOT/'v32_fallback_correction.py')==r['correction_code_sha256']
 elif version=='33':
  r=json.loads((old/'fallback_correction_amendment.json').read_text());assert m.sha(m.ROOT/'run_v33_fixed.py')==r['adapter_sha256'];assert m.sha(m.ROOT/'v32_fallback_correction.py')==r['correction_sha256']
 else:
  r=json.loads((old/'integration_correction_amendment.json').read_text());assert m.sha(m.ROOT/'v34_corrected_methods.py')==r['corrected_source_sha256'];r=json.loads((old/'count_projection_recovery.json').read_text());assert m.sha(m.ROOT/'run_v34_projected_counts.py')==r['recovery_code_sha256']
 with threadpool_limits(limits=2):
  if version=='32':
   import run_v32_experiment as run
   import v32_fallback_correction as correction
   m.A=dest;run.main();correction.main()
  elif version=='33':
   import run_v33_fixed as run
   import v32_fallback_correction as correction
   correction.install();run.A=dest;run.main()
  else:
   import v34_corrected_methods as f
   import run_v34_experiment as run
   import run_v34_projected_counts as counts
   f.A=dest;run.f=f;counts.f=f
   for outer in m.OUTERS:f.train(outer)
   run.mean_screen();counts.count_test('null');counts.count_test('treatment')
 print('Reproduction complete:',dest,flush=True)