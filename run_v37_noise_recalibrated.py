"""Final separately frozen combination; stop after nulls if its control gate fails."""
import argparse,json
from threadpoolctl import threadpool_limits
import v35_methods as m
import run_v37_recalibrated as x
from v37_noise_recalibrated import NAME,emit
from run_v37_ensemble import read,write
A=x.PARENT/'noise_recalibrated';x.A=A;x.NAME=NAME;x.emit=emit

def freeze():
 x.freeze();p=read(A/'plan.json');p['dependencies'].update({n:m.sha(m.ROOT/n) for n in ['v37_noise_recalibrated.py','run_v37_noise_recalibrated.py']})
 p.update({'arm':NAME,'known_outcomes':'All five primary and the post-mixture calibration development/null outcomes known. Noise mix passed all treatment guards but failed nulls; post-mixture calibration improved nulls but failed treatment guards. No outcome for this seventh fixed arm known.','formula':'Use the fixed mean50 response and initial0.30expected+0.65Poisson+0.05empirical population, then one further unchanged dual-moment calibration, same seed and200 iterations. No parent row mixing, weight grid or new fitted parameter.','early_stop':'Require every corrected control-null DE-call count <=V24. If any fails, do not run treatment counts or additional confirmation.'});write(A/'plan.json',p)

def gate():
 x.check();r=read(A/'null_results.json');base={v['family']:v['null_DE_calls'] for v in r['records'] if v['arm']=='v24'};rows=[{'family':v['family'],'v24':base[v['family']],'candidate':v['null_DE_calls'],'passes':v['null_DE_calls']<=base[v['family']]} for v in r['records'] if v['arm']==NAME];passed=all(v['passes'] for v in rows);write(A/'null_gate.json',{'selected':'pending' if passed else 'v24','passed':passed,'eligible':[],'gain_scores':{},'checks':{NAME:{'null_no_worse':passed}},'rows':rows,'null_results_sha256':m.sha(A/'null_results.json'),'treatment_counts_opened':False});print(json.dumps(read(A/'null_gate.json'),indent=2),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','null','gate','development','select','confirmation','confirm','audit']);a=p.parse_args()
 with threadpool_limits(limits=2):
  if a.stage=='freeze':freeze()
  elif a.stage=='null':x.null()
  elif a.stage=='gate':gate()
  elif a.stage=='development':assert read(A/'null_gate.json')['passed'];x.development()
  elif a.stage in ('select','confirm'):x.install();x.r.decide('development' if a.stage=='select' else 'confirmation')
  elif a.stage=='confirmation':x.install();x.r.count('confirmation')
  else:
   import audit_v37_ensemble as audit
   audit.A=A;audit.check=x.check;audit.main()
