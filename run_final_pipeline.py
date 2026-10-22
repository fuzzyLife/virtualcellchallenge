"""Execute final preparation in dependency order; never upload or push."""
from pathlib import Path
import time,sys,subprocess,json
R=Path(__file__).resolve().parent;A=R/'experiments/final_v91/analysis';PY=R/'.competition_env/Scripts/python.exe'
REQUIRED=['cd4_collection.json','H1_new_source.json','HCT116_new_source.json','HEK293T_new_source.json','validation_replay_A.json','validation_replay_B.json','validation_replay_C.json','count_runtime_replay.json','cli_final_smoke.json']
def run(script,*args):
 print('FINAL STAGE',script,*args,flush=True);subprocess.run([str(PY),'-X','utf8',str(R/script),*map(str,args)],cwd=R,check=True)
def main():
 start=time.monotonic();last=None
 while True:
  missing=[n for n in REQUIRED if not (A/n).exists()]
  if not missing:break
  if missing!=last:print('Final pipeline waiting for',', '.join(missing),flush=True);last=missing
  if time.monotonic()-start>6*3600:raise TimeoutError('Input collection/replay prerequisites did not finish. No upload attempted.')
  time.sleep(5)
 for name in ['validation_replay_A.json','validation_replay_B.json','validation_replay_C.json','count_runtime_replay.json','cli_final_smoke.json']:assert json.loads((A/name).read_text())['passed']
 run('build_final_sources.py');run('audit_final_sources.py');run('fit_final_submission.py','freeze');run('fit_final_submission.py','prior');run('fit_final_submission.py','prior','--shuffled');run('fit_final_submission.py','fit');run('audit_final_model.py');run('build_final_submission.py','--preflight');run('build_final_submission.py','--workers',4);run('package_final_submission.py');run('validate_final_submission.py');(A/'experiment_complete.json').write_text(json.dumps(dict(completed=True,partition='test',final_upload_ready=True,uploaded=False,final_accuracy_unknown=True),indent=2),encoding='utf-8');print('FINAL D/E/F SUBMISSION READY; NO UPLOAD ATTEMPTED',flush=True)
if __name__=='__main__':main()
