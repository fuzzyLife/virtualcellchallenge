"""Execute all prespecified H1 ideas, then expand only a passing locked finalist."""
import sys,json,subprocess
from pathlib import Path
from v27_methods import ROOT,A,load_plan,write
import build_submission as b

def call(name,*args):
 subprocess.run([sys.executable,str(ROOT/name),*args],cwd=ROOT,check=True)

def main():
 load_plan()
 if not (A/'sampling_evaluation.json').exists():call('run_v27_experiment.py','ablation')
 # Reproduce source sampling in the new experiment; leave old training arrays intact.
 call('train_v27_decoder.py','prepare','--family','H1')
 from train_v27_decoder import train
 train('H1','panel')
 completed=[]
 for key in ['H1','HCT116','HEK293T','CD4']:
  if key!='H1':call('train_v27_decoder.py','train','--family',key)
  for stage in ['null','treatment']:
   if not (A/f'{stage}_{key}_evaluation.json').exists():call('run_v27_experiment.py',stage,'--family',key)
   if not (A/f'{stage}_{key}_result.json').exists():call('assess_v27_experiment.py',stage,'--family',key)
   # H1 treatment is diagnostic even after null failure. For later families,
   # calibration failure still permits that prespecified paired diagnostic.
  completed.append(key)
  results={s:json.loads((A/f'{s}_{key}_result.json').read_text()) for s in ['null','treatment']}
  if not all(r['passed'] for r in results.values()):
   write(A/'release_decision.json',{'released':False,'uploaded':False,'new_package_generated':False,'reserve_opened':False,'completed_families':completed,'failed_family':key,'failed_checks':{s:[k for k,v in r['checks'].items() if not v] for s,r in results.items()},'reason':'Locked finalist failed prespecified control/treatment checks. All nine pilot arms were intervention-tested; no adaptive substitution. Preserve V24.','result_hashes':{s:b.sha256(A/f'{s}_{key}_result.json') for s in results},'plan_sha256':b.sha256(A/'plan.json')});print('V27 rejected after all authorized pilot comparisons',key,flush=True);return
 # A complete passing development run still needs a locked fresh-reserve test;
 # no generator can advertise an upload-ready artifact before that and audits.
 write(A/'development_pass.json',{'passed':True,'completed_families':completed,'plan_sha256':b.sha256(A/'plan.json'),'next_required':'Lock finalist, complete author-MD5-verified D3Rest reserve with ALL CD4 excluded, run guards, official prep and full-payload audits; no upload yet.'})
 print('V27 full development passed; reserve and release audit now required',flush=True)
if __name__=='__main__':main()
