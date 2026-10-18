"""Resume independent V57 context inferences with at most two GPU workers.

Frozen inference, cohorts and seeds are unchanged; distinct contexts write distinct files.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess,sys
from v57_pie import *

def run(script,*args):
 print('V57 stage:',script,*args,flush=True)
 subprocess.run([sys.executable,'-X','utf8',str(ROOT/script),*args],check=True)

def complete(source):
 p=A/'predictions'/(source+'.json');done=A/'predictions'/(source+'_evidence_ablation.json')
 if not p.exists() or not done.exists():return False
 r=read(p);path=CACHE/'predictions'/(source+'.npz')
 return r['plan_sha256']==sha(A/'plan.json') and all(str(i) in r['records'] for i in r['required_targets']) and r['sha256']==sha(path)

def one(source):
 if complete(source):print('Verified complete inference:',source,flush=True);return
 run('infer_v57_pie.py','--source',source)
 assert complete(source)

def main():
 assert read(A/'evidence_independent_audit.json')['passed'];assert all((A/'inputs'/(s+'.json')).exists() for s in SOURCES)
 with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,SOURCES))
 run('run_v57_cached_stage.py','fit_v57_pie_decoder.py','train');run('run_v57_cached_stage.py','fit_v57_pie_decoder.py','evaluate')
 write(A/'response_pipeline_complete.json',{'passed':True,'plan_sha256':sha(A/'plan.json'),'scope':'Response evaluation complete; count stage separate','orchestrator_sha256':sha(Path(__file__)),'maximum_concurrent_GPU_contexts':2,'scientific_settings_changed':False})
 print('Response stage complete. Inspect response_summary.csv before count generation.',flush=True)

if __name__=='__main__':main()
