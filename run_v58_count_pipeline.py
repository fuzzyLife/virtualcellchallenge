"""Run the frozen V58 count comparisons on independent contexts, then audit/report."""
from concurrent.futures import ThreadPoolExecutor
import subprocess,sys
from v58_pie import *

def run(script,*args):
 print('V58 count stage:',script,*args,flush=True)
 subprocess.run([sys.executable,'-X','utf8',str(ROOT/script),*args],check=True)

def main():
 assert read(A/'response_evaluation.json')['complete'], 'Inspect response results before invoking the count stage'
 with ThreadPoolExecutor(max_workers=2) as pool:
  list(pool.map(lambda s:run('run_v58_cached_stage.py','evaluate_v58_pie_counts.py','run','--source',s),SOURCES))
 run('evaluate_v58_pie_counts.py','summarize')
 run('run_v58_cached_stage.py','summarize_v58_panels.py')
 run('audit_v58_pie.py','artifacts')
 run('run_v58_cached_stage.py','document_v58_pie.py')
 write(A/'count_pipeline_complete.json',{'passed':True,'plan_sha256':sha(A/'plan.json'),'orchestrator_sha256':sha(ROOT/'run_v58_count_pipeline.py'),'maximum_concurrent_CPU_contexts':2,'scientific_settings_changed':False})

if __name__=='__main__':main()
