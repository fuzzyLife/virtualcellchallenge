"""Resume the frozen V57 pipeline through response evaluation, before count generation."""
from pathlib import Path
import subprocess,sys,time
from v57_pie import A,SOURCES,sha,read,write

def wait_for(path,seconds=3600):
 started=time.monotonic();last=-31
 while not path.exists():
  elapsed=time.monotonic()-started
  if elapsed>seconds:raise TimeoutError('Missing prerequisite: '+str(path))
  if elapsed-last>=30:print('Waiting for prerequisite:',path.name,flush=True);last=elapsed
  time.sleep(5)

def run(script,*args):
 print('V57 stage:',script,*args,flush=True)
 subprocess.run([sys.executable,'-X','utf8',str(Path(__file__).resolve().parent/script),*args],check=True)

def main():
 wait_for(A/'evidence_receipt.json')
 run('audit_v57_pie_cached.py','evidence')
 for source in SOURCES:wait_for(A/'inputs'/(source+'.json'))
 for source in SOURCES:run('infer_v57_pie.py','--source',source)
 run('fit_v57_pie_decoder.py','train')
 run('fit_v57_pie_decoder.py','evaluate')
 write(A/'response_pipeline_complete.json',{'passed':True,'plan_sha256':sha(A/'plan.json'),'scope':'Response evaluation complete; count stage intentionally separate','orchestrator_sha256':sha(Path(__file__))})
 print('Response stage complete. Inspect response_summary.csv before count generation.',flush=True)

if __name__=='__main__':main()
