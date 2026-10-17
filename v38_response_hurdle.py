"""V38: separate activity/direction, conditional magnitude, and count variation.

Missing response measurements are NaN features, never zero-response votes.
Whole biological studies are excluded from recipient features and outer fits.
"""
import gc, hashlib, json
from functools import lru_cache
import joblib
import numpy as np
from scipy.stats import rankdata
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.metrics import average_precision_score, log_loss
import build_submission as b
from fit_v7_response import read_cache
from v15_ensemble_methods import cube, pool, profiles
from v7_response_methods import FAMILIES

ROOT=b.ROOT
OUT=ROOT/'experiments/v38_response_hurdle'
A=OUT/'analysis'
SEED=38038
OUTERS=('H1','Orion','CD4')
ARMS=('hurdle25','weighted25','selective50')
GROUP={'K562':'Replogle','CD4':'CD4','H1':'H1','HCT116':'Orion','HEK293T':'Orion'}
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def write(p,v):b.write_json(p,v)
def split(t):return int(hashlib.sha256(('v38-target:'+str(t)).encode()).hexdigest()[:8],16)%5
def legacy(ex):return [str(f) for f in FAMILIES if GROUP[str(f)] in ex]
@lru_cache(1)
def data():return read_cache()
@lru_cache(1)
def embeddings():
 from v30_methods import embeddings as load
 names,values=load();mp={str(n):v[:8] for n,v in zip(names,values)}
 return np.array([mp.get(str(n),np.zeros(8)) for n in data()['genes']],np.float32)
@lru_cache(16)
def bank(excluded):
 d=data();y,v,ok,own=cube(d,legacy(excluded));pp=profiles(d)
 pooled=pool(y,v,ok,np.array([2.,.25,2.,1.,1.]),own)
 return y,v,ok,pp,own,pooled

def descriptor(query,excluded):
 y,v,ok,pp,own,pooled=bank(tuple(sorted(excluded)))
 q=np.asarray(query,float);qr=np.full(len(q),np.nan);qm=np.isfinite(q)
 qr[qm]=rankdata(q[qm])/qm.sum();ranked=[]
 for group in ('Replogle','CD4','H1','Orion'):
  ids=[i for i,f in enumerate(FAMILIES) if GROUP[str(f)]==group and group not in excluded]
  if not ids:continue
  vv=ok[ids];n=vv.sum(0);mean=np.divide(np.where(vv,y[ids],0).sum(0),n,out=np.zeros_like(pooled),where=n>0)
  var=np.divide(np.where(vv,v[ids]+(y[ids]-mean)**2,0).sum(0),n,out=np.ones_like(pooled)*1e6,where=n>0)
  available=n>0;pm=np.isfinite(pp[ids]);den=pm.sum(0)
  fp=np.divide(np.where(pm,pp[ids],0).sum(0),den,out=np.full(len(q),np.nan),where=den>0)
  good=np.isfinite(q)&np.isfinite(fp)&~np.isin(data()['genes'],data()['targets'])
  rho=float(np.corrcoef(rankdata(q[good]),rankdata(fp[good]))[0,1]) if good.sum()>100 else -1.
  ranked.append((rho,group,mean,np.sqrt(var),available,fp))
 ranked.sort(key=lambda x:(-x[0],x[1]))
 return {'q':q,'rank':qr,'sources':ranked,'pooled':pooled,'own':own,'excluded':list(excluded)}

def features(desc,ti,positions):
 pos=np.asarray(positions);q=desc['q'];o=desc['own'][ti];e=embeddings();eb=e[pos];et=e[o]
 z=[np.log1p(q[pos]),desc['rank'][pos],np.full(len(pos),np.log1p(q[o])),
    np.full(len(pos),desc['rank'][o]),np.sum(eb*et,axis=1),desc['pooled'][ti,pos]]
 sign=[];errors=[]
 for k in range(4):
  if k>=len(desc['sources']):
   z.extend([np.full(len(pos),np.nan) for _ in range(8)]);continue
  rho,group,y,se,ok,fp=desc['sources'][k];mask=ok[ti,pos]
  val=np.where(mask,y[ti,pos],np.nan);err=np.where(mask,se[ti,pos],np.nan)
  z.extend([val,abs(val),err,np.clip(val/err,-20,20),mask.astype(float),
            np.log1p(q[pos])-np.log1p(fp[pos]),np.full(len(pos),rho),
            np.full(len(pos),np.log1p(q[o])-np.log1p(fp[o]))])
  sign.append(np.where(mask,np.sign(y[ti,pos]),0));errors.append(mask)
 n=np.sum(errors,axis=0) if errors else np.zeros(len(pos))
 z.extend([n,np.divide(abs(np.sum(sign,axis=0)),n,out=np.zeros(len(pos)),where=n>0)])
 x=np.column_stack(z+[eb[:,i] for i in range(8)]+[np.full(len(pos),v) for v in et]).astype(np.float32)
 x[~np.isfinite(x)]=np.nan
 return x,n

def table(outer):
 d=data();xs=[];ys=[];labs=[];weights=[];tids=[];rids=[];prior=[]
 ex=[] if outer=='all' else [outer]
 counts={g:sum(GROUP[str(f)]==g for f in d['row_families']) for g in set(GROUP.values())}
 for ci in range(7):
  rows=np.flatnonzero(d['row_contexts']==ci)
  if not len(rows):continue
  group=GROUP[str(d['row_families'][rows[0]])]
  if group in ex:continue
  ds=descriptor(d['profiles'][ci],tuple(sorted(set(ex+[group]))))
  for ri in rows:
   ti=int(d['row_targets'][ri]);valid=d['valid'][ri]&np.isfinite(d['profiles'][ci])&(d['profiles'][ci]>5)
   valid[ds['own'][ti]]=False;possible=np.flatnonzero(valid)
   rng=np.random.default_rng(SEED+int(ri));pos=rng.choice(possible,min(256,len(possible)),replace=False)
   x,n=features(ds,ti,pos);keep=n>0;pos=pos[keep];x=x[keep]
   truth=d['truth'][ri,pos];se=d['errors'][ri,pos]
   active=(abs(truth)>=.25)&(abs(truth)>=2*se)
   label=np.where(active,np.where(truth>0,2,0),1).astype(np.int32)
   xs.append(x);ys.append(truth);labs.append(label);prior.append(ds['pooled'][ti,pos])
   weights.append(np.full(len(pos),1/max(counts[group],1)))
   tids.extend([ti]*len(pos));rids.extend([int(ri)]*len(pos))
  print('V38 fitting table',outer,ci,sum(len(x) for x in xs),flush=True)
 x=np.concatenate(xs);y=np.concatenate(ys);label=np.concatenate(labs);w=np.concatenate(weights);w/=w.mean()
 return x,y,label,w,np.asarray(tids),np.asarray(rids),np.concatenate(prior)

def fit(outer):
 folder=A/'models';folder.mkdir(parents=True,exist_ok=True);path=folder/(outer+'.joblib');assert not path.exists()
 x,y,label,w,ti,ri,prior=table(outer);val=np.array([split(data()['targets'][i])==0 for i in ti]);train=~val
 params=dict(max_iter=120,max_leaf_nodes=15,min_samples_leaf=80,l2_regularization=10.,learning_rate=.06,early_stopping=False,random_state=SEED)
 clf=HistGradientBoostingClassifier(**params).fit(x[train],label[train],sample_weight=w[train]);mags={}
 for lab in (0,2):
  ix=train&(label==lab);assert ix.sum()>=50,(outer,lab,int(ix.sum()))
  mags[lab]=HistGradientBoostingRegressor(**params).fit(x[ix],np.log1p(abs(y[ix])),sample_weight=w[ix])
 reg=HistGradientBoostingRegressor(**params).fit(x[train],y[train]-prior[train],sample_weight=w[train]*np.where(label[train]!=1,8.,1.))
 model={'classifier':clf,'magnitude':mags,'weighted':reg,'outer':outer,'class_rule':'abs(LFC)>=.25 and abs(LFC)>=2*SE; diagnostic only','seed':SEED}
 pp=clf.predict_proba(x[val]);hm=pp[:,2]*np.maximum(np.expm1(mags[2].predict(x[val])),0)-pp[:,0]*np.maximum(np.expm1(mags[0].predict(x[val])),0)
 wr=prior[val]+reg.predict(x[val]);vv=w[val];yy=y[val]
 diagnostics={'outer':outer,'training_pairs':int(train.sum()),'validation_pairs':int(val.sum()),'training_classes':np.bincount(label[train],minlength=3).tolist(),'validation_classes':np.bincount(label[val],minlength=3).tolist(),
  'activity_average_precision':float(average_precision_score(label[val]!=1,1-pp[:,1],sample_weight=vv)),
  'class_log_loss':float(log_loss(label[val],pp,sample_weight=vv,labels=[0,1,2])),
  'validation_mse':{n:float(np.average((p-yy)**2,weights=vv)) for n,p in [('pooled',prior[val]),('hurdle',hm),('weighted',wr)]},
  'validation_active_mae':{n:float(np.average(abs(p[label[val]!=1]-yy[label[val]!=1]),weights=vv[label[val]!=1])) for n,p in [('pooled',prior[val]),('hurdle',hm),('weighted',wr)]},
  'recipient_studies':[g for g in sorted(set(GROUP.values())) if outer=='all' or g!=outer],
  'feature_exclusions':'outer entire study plus each training recipient entire study','validation_split':'All occurrences of hash0 targets excluded from fitting; allowed experimental donor outcomes for those targets remain retrieval features.'}
 joblib.dump(model,path,compress=3);write(path.with_suffix('.json'),{**diagnostics,'sha256':b.sha256(path),'plan_sha256':b.sha256(A/'plan.json')});print(json.dumps(diagnostics,indent=2),flush=True)
 del x,model;gc.collect()

def load(outer):
 p=A/'models'/(outer+'.joblib');rec=read(p.with_suffix('.json'));assert b.sha256(p)==rec['sha256'];assert rec['plan_sha256']==b.sha256(A/'plan.json')
 return joblib.load(p)

def predict(model,query,base,targets=None):
 d=data();ts=d['targets'] if targets is None else np.asarray(targets);ex=() if model['outer']=='all' else (model['outer'],);desc=descriptor(query,ex);pos=np.flatnonzero(np.isfinite(query)&(query>1));result={a:base.copy() for a in ARMS};trace=[]
 for j,t in enumerate(ts):
  ti=int(np.flatnonzero(d['targets']==t)[0]);x,n=features(desc,ti,pos);pp=model['classifier'].predict_proba(x)
  down=np.clip(np.expm1(model['magnitude'][0].predict(x)),0,3);up=np.clip(np.expm1(model['magnitude'][2].predict(x)),0,3)
  hm=pp[:,2]*up-pp[:,0]*down;delta=np.clip(hm-desc['pooled'][ti,pos],-1,1)
  rd=np.clip(model['weighted'].predict(x),-1,1)
  active=1-pp[:,1];agreement=np.divide(abs(pp[:,2]-pp[:,0]),active,out=np.zeros_like(active),where=active>0)
  gate=(active>=.20)&(agreement>=.8)&(hm*desc['pooled'][ti,pos]>0)&(abs(hm)>abs(desc['pooled'][ti,pos]))
  result['hurdle25'][j,pos]+=.25*np.where(n>0,delta,0)
  result['weighted25'][j,pos]+=.25*np.where(n>0,rd,0)
  result['selective50'][j,pos]+=.5*np.where((n>0)&gate,delta,0)
  trace.append({'target':str(t),'source_supported_outputs':int((n>0).sum()),'selective_outputs':int(gate.sum()),'maximum_activity':float(active.max(initial=0))})
 for a in ARMS:
  result[a]=np.clip(result[a],-3,3).astype(np.float32)
  for j,t in enumerate(ts):result[a][j,d['genes']==t]=0
 return result,trace

def emit(template,effect,own,seed,prepared,atlas,transport=False):
 from v25_methods import emit as standard
 count,diag=standard(template,effect,own,seed,prepared,atlas)
 if transport:
  from v27_methods import anchored_initial
  from v29_methods import project
  count,extra=project(anchored_initial(template,effect,own,prepared[1]),count,seed+270270);diag.update(extra)
 return count,diag

