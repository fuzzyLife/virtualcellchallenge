"""Run V58 corrected source-only fitting and response evaluation before counts."""
import subprocess,sys
from v58_pie import *

def run(script,*args):
 print('V58 stage:',script,*args,flush=True);subprocess.run([sys.executable,'-X','utf8',str(ROOT/script),*args],check=True)

def main():
 run('prepare_v58_pie.py','prepare');run('audit_v58_pie.py','boundaries')
 run('run_v58_cached_stage.py','fit_v58_pie_decoder.py','train')
 run('run_v58_cached_stage.py','fit_v58_pie_decoder.py','evaluate')
 write(A/'response_pipeline_complete.json',{'passed':True,'plan_sha256':sha(A/'plan.json'),'orchestrator_sha256':sha(Path(__file__)),'count_stage_started':False})

if __name__=='__main__':main()
