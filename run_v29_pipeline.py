"""Gated V29 follow-up, preserving all prior numeric experiments."""
import sys,json,subprocess
from v29_methods import ROOT,A,load_plan,write
import build_submission as b

def call(name,*args):subprocess.run([sys.executable,str(ROOT/name),*args],cwd=ROOT,check=True)
def reject(stage,key,completed):
 r=json.loads((A/f'{stage}_{key}_result.json').read_text());write(A/'release_decision.json',{'released':False,'uploaded':False,'new_package_generated':False,'reserve_opened':False,'failed_stage':stage,'failed_family':key,'completed':completed,'failed_checks':[k for k,v in r['checks'].items() if not v],'result_sha256':b.sha256(A/f'{stage}_{key}_result.json'),'plan_sha256':b.sha256(A/'plan.json'),'reason':'Sole global-transport finalist failed frozen checks; no diagnostic substitution or reserve opening. Preserve V24.'});print('V29 rejected',stage,key,flush=True)
def main():
 load_plan();completed={'null':[],'treatment':[]}
 for stage in ['null','treatment']:
  for key in ['H1','HCT116','HEK293T','CD4']:
   if not (A/f'{stage}_{key}_evaluation.json').exists():call('run_v29_experiment.py',stage,'--family',key)
   if not (A/f'{stage}_{key}_result.json').exists():call('assess_v29_experiment.py',stage,'--family',key)
   completed[stage].append(key)
   if not json.loads((A/f'{stage}_{key}_result.json').read_text())['passed']:reject(stage,key,completed);return
 if not (A/'development_result.json').exists():call('assess_v29_experiment.py','development')
 r=json.loads((A/'development_result.json').read_text())
 if not r['passed']:
  write(A/'release_decision.json',{'released':False,'uploaded':False,'new_package_generated':False,'reserve_opened':False,'failed_stage':'development','completed':completed,'failed_checks':[k for k,v in r['checks'].items() if not v],'result_sha256':b.sha256(A/'development_result.json'),'plan_sha256':b.sha256(A/'plan.json'),'reason':'Global development criteria failed; no reserve opening or new package. Preserve V24.'});return
 write(A/'development_pass.json',{'passed':True,'plan_sha256':b.sha256(A/'plan.json'),'result_sha256':b.sha256(A/'development_result.json'),'next_required':'Lock finalist, complete MD5-verified fresh D3Rest reserve excluding all CD4 response information, then official packaging and full-payload audits.'});print('V29 development passed; fresh reserve and audits required',flush=True)
if __name__=='__main__':main()
