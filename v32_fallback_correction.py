"""V32 correction: retain exact V24 fallback outside supported response tokens."""
import json
import numpy as np
import v32_methods as m
from run_v32_experiment import evaluate
from threadpoolctl import threadpool_limits

def predict(g,query,alpha,parent=None):
 data=m.old.read_cache();base=m.old.baseline(data,query,m.old.excluded_legacy(g['excluded'])) if parent is None else parent.copy();est,mask,routes=m.estimate(g,query);own=np.array([np.flatnonzero(data['genes']==t)[0] for t in data['targets']]);full=np.zeros_like(base);full[:,g['axis']]=est;full=np.array([m.old.gated(v.copy(),np.nan_to_num(query),int(o)) for v,o in zip(full,own)]);delta=np.zeros_like(base,dtype=float);delta[:,g['axis']]=np.where(mask,np.clip(full[:,g['axis']]-base[:,g['axis']],-1,1),0);pred=np.clip(base+alpha*delta,-3,3);pred[np.arange(300),own]=0;return pred,base,routes

def install():m.predict=predict

def load(path):
 rec=json.loads(path.with_suffix('.json').read_text());assert m.sha(path)==rec['sha256']
 with np.load(path,allow_pickle=False) as z:g={k:z[k].copy() for k in z.files}
 g.update({k:v for k,v in rec.items() if k not in ['sha256','plan_sha256']});keys,e=m.old.embeddings();g['embedding']={str(k):v[:16] for k,v in zip(keys,e)};return g

def main():
 import v32_runtime as runtime
 m.check();runtime.install();install();assert not (m.A/'corrected_results.json').exists();records=[];choices={}
 for outer in m.OUTERS:
  vals=[]
  for inner in m.OUTERS:
   if inner==outer:continue
   g=load(m.A/'models'/f'{outer}_{inner}.npz')
   for alpha in m.ALPHAS:
    r=evaluate(g,inner,alpha);records.append({'stage':'inner','outer':outer,'test':inner,'arm':'paired','alpha':alpha,**r});vals.append((alpha,r['mse_ratio']))
  alpha=min(m.ALPHAS,key=lambda a:(np.mean([v for aa,v in vals if aa==a]),a));choices[outer]=alpha
  for arm,suffix in [('paired',''),('identity_shuffled','_shuffle'),('source_copy','')]:
   g=load(m.A/'models'/f'{outer}{suffix}.npz')
   if arm=='source_copy':g['coef']=np.zeros_like(g['coef'])
   r=evaluate(g,outer,alpha);records.append({'stage':'outer','outer':outer,'test':outer,'arm':arm,'alpha':alpha,**r});print('V32 corrected',outer,arm,r['mse_ratio'],r['mae_ratio'],flush=True)
 summary={arm:{k:float(np.mean([r[k] for r in records if r['stage']=='outer' and r['arm']==arm])) for k in ['mse_ratio','mae_ratio']} for arm in ['paired','identity_shuffled','source_copy']};s=summary['paired'];checks={'MSE_gain':s['mse_ratio']<=.99,'MAE_gain':s['mae_ratio']<=.995,'worst_family':max(r['mse_ratio'] for r in records if r['stage']=='outer' and r['arm']=='paired')<=1.02,'beats_shuffle':all(s[k]<summary['identity_shuffled'][k] for k in s),'beats_copy':all(s[k]<summary['source_copy'][k] for k in s)};m.write(m.A/'corrected_results.json',{'parent_plan_sha256':m.sha(m.A/'plan.json'),'original_results_sha256':m.sha(m.A/'results.json'),'correction_code_sha256':m.sha(m.ROOT/'v32_fallback_correction.py'),'reason':'Original code reapplied query-expression gating to already-gated V24 fallback, including outside modeled outputs. Correction changes only supported source-token coordinates. Original artifacts preserved; saved fits replayed, same inner alpha grid and thresholds. No new model/hyperparameter search.','records':records,'choices':choices,'summary':summary,'checks':checks,'passes':all(checks.values())});print(json.dumps({'summary':summary,'checks':checks},indent=2))
if __name__=='__main__':
 with threadpool_limits(limits=2):main()