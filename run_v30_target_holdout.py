"""Separately frozen within-source unseen-target stress test; never a release selector."""
from pathlib import Path
import json,hashlib
import numpy as np
from sklearn.utils.extmath import randomized_svd
from threadpoolctl import threadpool_limits
import v30_methods as m
import build_submission as b
D=m.A/'target_holdout'

def fit(source,ids,variant):
 d=source;keys,e=m.embeddings();lookup={t:i for i,t in enumerate(keys)}
 if variant=='random':e=e[np.random.default_rng(30301).permutation(len(e))]
 names=d['targets'][ids];x=e[[lookup[t] for t in names]];mu=x.mean(0);sd=np.maximum(x.std(0),.01);x=(x-mu)/sd;y=d['y'][ids].astype(float);w=d['n'][ids].astype(float);w=w/(w+100);w/=w.mean();gene_names=m.read_cache()['genes'][d['positions']];gp={g:i for i,g in enumerate(gene_names)};own=np.array([gp.get(t,-1) for t in names]);mask=np.ones(y.shape,bool);ok=np.flatnonzero(own>=0);mask[ok,own[ok]]=False;y[~mask]=0;bias=(w[:,None]*y).sum(0)/(w[:,None]*mask).sum(0);center=y-bias;center[~mask]=0;scale=np.maximum(np.sqrt((w[:,None]*center**2).sum(0)/(w[:,None]*mask).sum(0)),.1);center/=scale;_,_,vt=randomized_svd(center*np.sqrt(w[:,None]),n_components=32,random_state=30301,n_iter=4);g=vt.T;models={}
 for r in [16,32]:
  for lam in [1.,10.]:models[(r,lam)]=m.masked_bilinear(x,g[:,:r],center,own,w,lam*len(ids)/64)
 return {'basis':g,'bias':bias,'scale':scale,'feature_mean':mu,'feature_scale':sd,'weights':models,'embedding':e,'lookup':lookup,'training_targets':names}

def score(d,ids,model,config):
 r,lam=config;names=d['targets'][ids];x=(model['embedding'][[model['lookup'][t] for t in names]]-model['feature_mean'])/model['feature_scale'];pred=np.clip((x@model['weights'][config]@model['basis'][:,:r].T)*model['scale']+model['bias'],-3,3);truth=d['y'][ids];mask=np.ones(truth.shape,bool);gp={g:i for i,g in enumerate(m.read_cache()['genes'][d['positions']])}
 for i,t in enumerate(names):
  if t in gp:mask[i,gp[t]]=False
 error=(pred-truth)[mask];constant=(np.broadcast_to(model['bias'],truth.shape)-truth)[mask];zero=truth[mask];return {'mse_vs_constant':float((error@error)/(constant@constant)),'mae_vs_constant':float(abs(error).sum()/abs(constant).sum()),'mse_vs_zero':float((error@error)/(zero@zero)),'targets':len(ids),'coordinates':int(mask.sum())}

def main():
 m.load_plan();D.mkdir(exist_ok=True);plan=D/'plan.json'
 assert not (D/'results.json').exists(),'Stress-test outcomes already exist; do not reopen for selection'
 specs={}
 keys,_=m.embeddings()
 for name in ['K562_full','RPE1_full']:
  d=np.load(m.A/'sources'/f'{name}.npz',allow_pickle=False);names=d['targets'];eligible=np.isin(names,keys);part=np.array([int(hashlib.sha256(('v30-target:'+str(t)).encode()).hexdigest()[:8],16)%100 for t in names]);specs[name]={'train':np.flatnonzero(eligible&(part<60)).tolist(),'inner':np.flatnonzero(eligible&(part>=60)&(part<80)).tolist(),'test':np.flatnonzero(eligible&(part>=80)).tolist(),'source_sha256':b.sha256(m.A/'sources'/f'{name}.npz')}
 body={'parent_plan_sha256':b.sha256(m.A/'plan.json'),'script_sha256':b.sha256(Path(__file__)),'sources':specs,'variants':['string','random'],'rank':[16,32],'lambda':[1,10],'selection':'Within-source 60/20/20 target-identity hash split, same target cannot enter multiple roles. Choose hyperparameters on inner MSE versus fitting-target mean; refit on80%, then score20% once. This diagnoses unseen-target learning and cannot rescue a failed context-transfer release candidate.','limits':'Source controls and source eligibility metadata used; all test treatment effects excluded from fitting/bases and model selection. This is within-source transfer, not an unseen-lineage test. It is independently frozen after V30 H1 context outcomes were visible; no target-stress outcome was opened before this plan.'}
 if plan.exists():assert json.loads(plan.read_text())==body
 else:m.write(plan,body)
 results=[]
 for name,spec in specs.items():
  d=np.load(m.A/'sources'/f'{name}.npz',allow_pickle=False)
  for variant in ['string','random']:
   first=fit(d,np.array(spec['train']),variant);inner={str(k):score(d,np.array(spec['inner']),first,k) for k in first['weights']};config=min(first['weights'],key=lambda k:(inner[str(k)]['mse_vs_constant'],k));del first;final=fit(d,np.array(spec['train']+spec['inner']),variant);test=score(d,np.array(spec['test']),final,config);rec={'source':name,'variant':variant,'config':list(config),'inner_results':inner,'test':test,'training_targets':final['training_targets'].tolist(),'test_target_names':d['targets'][spec['test']].tolist()};results.append(rec);file=D/f'{name}_{variant}.npz';np.savez_compressed(file,basis=final['basis'],bias=final['bias'],scale=final['scale'],feature_mean=final['feature_mean'],feature_scale=final['feature_scale'],weights=final['weights'][config]);rec['model_sha256']=b.sha256(file);print('V30 unseen-target stress',name,variant,test,flush=True)
 m.write(D/'results.json',{'plan_sha256':b.sha256(plan),'records':results,'release_eligible':False})
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
