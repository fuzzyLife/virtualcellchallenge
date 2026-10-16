"""V33 TabICLv2 versus matched CatBoost/ridge on paired response programs.

Local adaptation of tabular-foundation-model prediction, not paper reproduction.
Query contexts supply controls only. In-context examples come from public sources.
"""
import os,sys,json,argparse,gc,time,hashlib
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8')
from pathlib import Path
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'.v32_dependencies'))
import numpy as np
from threadpoolctl import threadpool_limits
from sklearn.linear_model import Ridge
import v32_methods as m
import v32_runtime as runtime
from run_v30_experiment import rows_for,summarize_context
A=ROOT/'experiments/v33_tabular_response/analysis';SEED=33033

def latent_labels(g,y,own):
 # Solve coefficients on observed coordinates; a missing own label is not zero.
 z=(y/g['scale'])@g['basis'];ok=np.flatnonzero(own>=0)
 if len(ok):
  v=g['basis'][own[ok]];z[ok]+=np.sum(z[ok]*v,1)[:,None]*v/np.maximum(1-np.sum(v*v,1),1e-8)[:,None]
 return z

def cap_rows(d):
 fit=np.flatnonzero(~d['validation']);rng=np.random.default_rng(SEED);routes=sorted(set(d['routes'][fit]));chunks=[list(rng.permutation(fit[d['routes'][fit]==r])) for r in routes];ids=[]
 for i in range(max(map(len,chunks))):
  for c in chunks:
   if i<len(c):ids.append(c[i])
 return np.array(ids[:512],int)

def freeze():
 assert not (A/'plan.json').exists();receipt=json.loads((ROOT/'reference/v33_tabular/receipt.json').read_text());m.write(A/'checkpoint_receipt.json',receipt);deps=['run_v33_experiment.py','test_v33_tabular.py','collect_v33_checkpoint.py','v32_runtime.py','v32_methods.py'];m.write(A/'plan.json',{'version':33,'dependencies':{n:m.sha(ROOT/n) for n in deps},'parent_plan_sha256':m.sha(m.A/'plan.json'),'checkpoint':receipt,'arms':['tabicl','catboost','ridge','identity_shuffled_tabicl'],'response_rank':8,'max_prompt_rows':512,'selection':'Fixed alpha0.5 primary; full alpha1 diagnostic only. No choosing a model or alpha on outer outcomes. Whole H1/Orion/CD4 exclusions; global target-hash20% withheld from all source fitting/basis labels. Same allowed source response tokens at query.','tabicl':{'n_estimators':2,'batch_size':1,'norm_methods':['none'],'device':'cuda','use_amp':False,'use_fa3':False,'kv_cache':True,'checkpoint_version':'tabicl-regressor-v2-20260212.ckpt'},'catboost':{'iterations':200,'depth':5,'learning_rate':.05,'loss':'MultiRMSE'},'ridge_alpha':10.,'labels':'Own-coordinate-excluded least-squares program coefficients; no own-target expression supplied. Source missing outputs restricted to shared measured axis, others fall back to V24.','gate':'Primary TabICL MSE<=0.99, MAE<=0.995 vs V24; worst-family MSE<=1.02; better than identity-shuffled TabICL on both. CatBoost and ridge are matched descriptive controls, cannot replace a failed primary after outer scoring.','scope':'Small-budget512-row8-program pilot of tabular foundation modeling, not exhaustive family tuning or published VCC score. Reused public biology. No reserve/package/upload.'})

def check():
 p=json.loads((A/'plan.json').read_text())
 for n,h in p['dependencies'].items():assert m.sha(ROOT/n)==h,n
 assert m.sha(ROOT/'reference/v33_tabular'/p['checkpoint']['url'].split('/')[-1])==p['checkpoint']['sha256'];m.check();return p

def main():
 from tabicl import TabICLRegressor
 from catboost import CatBoostRegressor
 import torch
 plan=check();runtime.install();assert not (A/'results.json').exists();(A/'predictions').mkdir(exist_ok=True);records=[];inventories=[]
 for outer in m.OUTERS:
  g=m.fit_geometry([outer]);g['basis']=g['basis'][:,:8];d=m.paired_table(g);g.update(feature_mean=d['feature_mean'],feature_scale=d['feature_scale']);ids=cap_rows(d);trainx=d['x'][ids].astype(np.float32);z=latent_labels(g,d['y'][ids],d['own'][ids]);ci_values=np.unique(m.old.read_cache()['row_contexts'][rows_for(m.old.read_cache(),outer)]);captured=[];original=m.predict_delta
  def capture(gg,x):captured.append(x.copy());return np.zeros((len(x),len(gg['axis'])))
  m.predict_delta=capture
  for ci in ci_values:m.estimate(g,m.old.read_cache()['profiles'][ci])
  m.predict_delta=original
  queryx=np.concatenate(captured).astype(np.float32);query_keys=[hashlib.sha256(v.tobytes()).hexdigest() for v in queryx];va=np.flatnonzero(d['validation']);allx=np.concatenate([d['x'][va].astype(np.float32),queryx]);nva=len(va);outputs={};component_log=[]
  print('V33 outer begin',outer,'fit rows',len(ids),'query tokens',len(queryx),'features',trainx.shape[1],flush=True)
  for arm in ['ridge','catboost','tabicl','identity_shuffled_tabicl']:
   yy=z.copy()
   if arm=='identity_shuffled_tabicl':
    rng=np.random.default_rng(SEED)
    for route in sorted(set(d['routes'][ids])):
     ii=np.flatnonzero(d['routes'][ids]==route);yy[ii]=yy[rng.permutation(ii)]
   if arm=='ridge':model=Ridge(alpha=10.);model.fit(trainx,yy);code=model.predict(allx)
   elif arm=='catboost':model=CatBoostRegressor(iterations=200,depth=5,learning_rate=.05,loss_function='MultiRMSE',random_seed=SEED,thread_count=2,verbose=False,allow_writing_files=False);model.fit(trainx,yy);code=model.predict(allx)
   else:
    model=TabICLRegressor(**plan['tabicl'],model_path=ROOT/'reference/v33_tabular/tabicl-regressor-v2-20260212.ckpt',allow_auto_download=False,random_state=SEED,n_jobs=2);code=np.zeros((len(allx),8))
    for j in range(8):
     start=time.time();model.fit(trainx,yy[:,j]);code[:,j]=model.predict(allx);component_log.append({'arm':arm,'component':j,'seconds':time.time()-start});print('V33',outer,arm,'component',j+1,'/8 seconds',round(time.time()-start,1),flush=True)
   delta=code@g['basis'].T*g['scale'];outputs[arm]=delta;lookup={k:v for k,v in zip(query_keys,delta[nva:])}
   def replay(gg,x):return np.array([lookup[hashlib.sha256(np.asarray(v,dtype=np.float32).tobytes()).hexdigest()] for v in x])
   m.predict_delta=replay
   for alpha in [.5,1.]:
    ctx=[]
    for ci in ci_values:
     rr=rows_for(m.old.read_cache(),outer);rr=rr[m.old.read_cache()['row_contexts'][rr]==ci];pred,base,routes=m.predict(g,m.old.read_cache()['profiles'][ci],alpha);ctx.append({**summarize_context(m.old.read_cache(),rr,pred,base),'context':int(ci),'sources':routes})
    record={'outer':outer,'arm':arm,'alpha':alpha,'mse_ratio':float(np.mean([r['mse_ratio'] for r in ctx])),'mae_ratio':float(np.mean([r['mae_ratio'] for r in ctx])),'contexts':ctx};records.append(record);print('V33 result',outer,arm,alpha,record['mse_ratio'],record['mae_ratio'],flush=True)
   m.predict_delta=original;del model;gc.collect();torch.cuda.empty_cache()
  file=A/'predictions'/f'{outer}.npz';np.savez_compressed(file,axis=g['axis'],basis=g['basis'],scale=g['scale'],trainx=trainx,trainy=z,fit_indices=ids,validation_indices=va,queryx=queryx,**outputs);inv={'outer':outer,'excluded_families':[outer],'fitting_sources':g['sources'],'fit_targets':d['targets'][ids].tolist(),'validation_targets':d['targets'][va].tolist(),'basis_training':g['basis_training'],'query_contexts':ci_values.tolist(),'sha256':m.sha(file),'components':component_log};m.write(file.with_suffix('.json'),inv);inventories.append(inv)
 summary={arm:{k:float(np.mean([r[k] for r in records if r['arm']==arm and r['alpha']==.5])) for k in ['mse_ratio','mae_ratio']} for arm in plan['arms']};s=summary['tabicl'];checks={'MSE_gain':s['mse_ratio']<=.99,'MAE_gain':s['mae_ratio']<=.995,'worst_family':max(r['mse_ratio'] for r in records if r['arm']=='tabicl' and r['alpha']==.5)<=1.02,'beats_shuffle':all(s[k]<summary['identity_shuffled_tabicl'][k] for k in s)};m.write(A/'results.json',{'plan_sha256':m.sha(A/'plan.json'),'records':records,'summary':summary,'checks':checks,'passes':all(checks.values()),'scope':'Public masked LFC mean ratios; fixed alpha0.5 primary, alpha1 diagnostic.','runtime':{'torch':torch.__version__,'GPU':torch.cuda.get_device_name(0)}});print(json.dumps({'summary':summary,'checks':checks},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','run']);a=p.parse_args()
 with threadpool_limits(limits=2):freeze() if a.stage=='freeze' else main()