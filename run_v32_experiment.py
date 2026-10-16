"""Nested held-family V32 response-transfer pilot."""
import json,argparse
import numpy as np
from threadpoolctl import threadpool_limits
import v32_methods as m
from run_v30_experiment import rows_for,summarize_context

def evaluate(g,key,alpha):
 d=m.old.read_cache();rows=rows_for(d,key);recs=[]
 for ci in np.unique(d['row_contexts'][rows]):
  rr=rows[d['row_contexts'][rows]==ci];p,b,routes=m.predict(g,d['profiles'][ci],alpha);recs.append({**summarize_context(d,rr,p,b),'context':int(ci),'sources':routes})
 return {k:float(np.mean([r[k] for r in recs])) for k in ['mse_ratio','mae_ratio']}|{'contexts':recs}
def main():
 m.check();assert not (m.A/'results.json').exists();(m.A/'models').mkdir(exist_ok=True);records=[];choices={}
 for outer in m.OUTERS:
  inners=[]
  for inner in m.OUTERS:
   if inner==outer:continue
   g=m.fit([outer,inner]);path=m.A/'models'/f'{outer}_{inner}.npz';m.save_model(g,path)
   for a in m.ALPHAS:
    r=evaluate(g,inner,a);records.append(dict(stage='inner',outer=outer,test=inner,arm='paired',alpha=a,**r));inners.append((a,r['mse_ratio']))
   print('V32 inner',outer,inner,flush=True)
  a=min(m.ALPHAS,key=lambda a:(np.mean([v for aa,v in inners if aa==a]),a));choices[outer]=a;g=m.fit([outer]);m.save_model(g,m.A/'models'/f'{outer}.npz');r=evaluate(g,outer,a);records.append(dict(stage='outer',outer=outer,test=outer,arm='paired',alpha=a,**r));print('V32 outer',outer,a,r['mse_ratio'],r['mae_ratio'],flush=True)
  sh=m.fit([outer],True);m.save_model(sh,m.A/'models'/f'{outer}_shuffle.npz');r=evaluate(sh,outer,a);records.append(dict(stage='outer',outer=outer,test=outer,arm='identity_shuffled',alpha=a,**r))
  cp=dict(g);cp['coef']=np.zeros_like(g['coef']);r=evaluate(cp,outer,a);records.append(dict(stage='outer',outer=outer,test=outer,arm='source_copy',alpha=a,**r))
 summary={arm:{k:float(np.mean([r[k] for r in records if r['stage']=='outer' and r['arm']==arm])) for k in ['mse_ratio','mae_ratio']} for arm in ['paired','identity_shuffled','source_copy']};s=summary['paired'];checks={'MSE_gain':s['mse_ratio']<=.99,'MAE_gain':s['mae_ratio']<=.995,'worst_family':max(r['mse_ratio'] for r in records if r['stage']=='outer' and r['arm']=='paired')<=1.02,'beats_shuffle':all(s[k]<summary['identity_shuffled'][k] for k in s),'beats_copy':all(s[k]<summary['source_copy'][k] for k in s)};m.write(m.A/'results.json',dict(plan_sha256=m.sha(m.A/'plan.json'),records=records,choices=choices,summary=summary,checks=checks,passes=all(checks.values()),scope='Masked public summary LFC error ratios versus held-family V24.'));print(json.dumps({'summary':summary,'checks':checks},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','run']);a=p.parse_args()
 with threadpool_limits(limits=2):m.freeze() if a.stage=='freeze' else main()