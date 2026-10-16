"""V32 paired-intervention transfer with control-conditioned shared response maps.

No recipient treatment is available at inference. Source-response tokens refer
only to the same target in allowed public sources. Never modifies V30 assets.
"""
from pathlib import Path
from functools import lru_cache
import hashlib,json,sys,time
import numpy as np
from scipy.stats import rankdata
from sklearn.utils.extmath import randomized_svd
from threadpoolctl import threadpool_limits
import v30_methods as old
import build_submission as b
ROOT=b.ROOT; A=ROOT/'experiments/v32_paired_context/analysis'; SEED=32032
OUTERS=['H1','Orion','CD4']; ALPHAS=[.25,.5,1.]

def write(p,x):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def sha(p):return b.sha256(p)
def family(s):return old.FAMILY[s]
def split(t):return int(hashlib.sha256(('v32-target:'+str(t)).encode()).hexdigest()[:8],16)%5
@lru_cache(None)
def sources():
 out={}
 for s in old.SOURCES:
  p=old.A/'sources'/f'{s}.npz';r=json.loads(p.with_suffix('.json').read_text());assert sha(p)==r['sha256']
  with np.load(p,allow_pickle=False) as z:out[s]={k:z[k].copy() for k in ['targets','positions','y','profile']}
 return out

def control_features(profile,axis):
 v=np.nan_to_num(profile[axis],nan=0.);return (rankdata(v,method='average')/len(v)).astype(np.float64)

def fit_geometry(excluded):
 data=old.read_cache();ss={s:d for s,d in sources().items() if family(s) not in excluded};assert len(ss)>=2
 common=set.intersection(*[set(d['positions']) for d in ss.values()]);common=np.array(sorted(common));assert len(common)>=512
 blocks=[];training=[]
 for s,d in ss.items():
  ids=np.flatnonzero([split(t)!=0 for t in d['targets']]);ids=np.random.default_rng(SEED).permutation(ids)[:512];lookup={int(v):i for i,v in enumerate(d['positions'])};p=[lookup[int(g)] for g in common];y=d['y'][ids][:,p].astype(float)
  ownmap={str(g):i for i,g in enumerate(data['genes'][common])}
  for i,t in enumerate(d['targets'][ids]):
   if str(t) in ownmap:y[i,ownmap[str(t)]]=0
  blocks.append(y);training.extend([s+':'+str(t) for t in d['targets'][ids]])
 # Equal source contribution in feature selection; no held-family labels.
 energy=np.mean([np.mean(v*v,0) for v in blocks],0);take=np.argsort(energy,kind='stable')[-min(2048,len(common)):];axis=common[take];y=np.concatenate([v[:,take]/np.sqrt(len(v)) for v in blocks]);scale=np.maximum(np.sqrt(np.mean([np.mean(v[:,take]**2,0) for v in blocks],0)),.1);_,_,vt=randomized_svd(y/scale,n_components=16,random_state=SEED,n_iter=4);basis=vt.T
 # Context coordinates are fit from allowed controls only, excluding panel genes.
 ca=common[~np.isin(data['genes'][common],data['targets'])];cf=np.array([control_features(d['profile'],ca) for d in ss.values()]);cm=cf.mean(0);_,_,cv=np.linalg.svd(cf-cm,full_matrices=False);cb=cv[:min(4,len(ss)-1)].T;cs=np.maximum(np.std((cf-cm)@cb,0),.01)
 names,vectors=old.embeddings();emb={str(k):v[:16] for k,v in zip(names,vectors)}
 return dict(sources=list(ss),excluded=list(excluded),axis=axis,scale=scale,basis=basis,context_axis=ca,context_mean=cm,context_basis=cb,context_scale=cs,embedding=emb,basis_training=training)

def ctx(g,profile):return np.clip((control_features(profile,g['context_axis'])-g['context_mean'])@g['context_basis']/g['context_scale'],-4,4)
def feature(g,e,source_profile,query_profile,target):
 code=(e/g['scale'])@g['basis'];s=ctx(g,source_profile);q=ctx(g,query_profile);delta=q-s;emb=g['embedding'].get(str(target),np.zeros(16))
 return np.r_[1.,code,s,q,np.outer(code,delta).ravel(),emb,np.outer(emb,delta).ravel()]

def paired_table(g,shuffle=False):
 ss=sources();data=old.read_cache();gx={str(v):i for i,v in enumerate(data['genes'][g['axis']])};xx=[];yy=[];ww=[];own=[];roles=[];routes=[];targets=[]
 for s in g['sources']:
  ds=ss[s];sm={str(t):i for i,t in enumerate(ds['targets'])};sp={int(v):i for i,v in enumerate(ds['positions'])};spi=[sp[int(v)] for v in g['axis']]
  for r in g['sources']:
   if s==r:continue
   dr=ss[r];rm={str(t):i for i,t in enumerate(dr['targets'])};rp={int(v):i for i,v in enumerate(dr['positions'])};rpi=[rp[int(v)] for v in g['axis']];shared=sorted(set(sm)&set(rm));train=[t for t in shared if split(t)!=0];valid=[t for t in shared if split(t)==0];rng=np.random.default_rng(SEED);train=list(np.array(train)[rng.permutation(len(train))[:256]]);valid=list(np.array(valid)[rng.permutation(len(valid))[:64]])
   if len(train)<8:continue
   order=train+valid;label_order=list(order)
   if shuffle:
    label_order=list(np.array(train)[rng.permutation(len(train))])+list(np.array(valid)[rng.permutation(len(valid))])
   for t,lt in zip(order,label_order):
    e=ds['y'][sm[t],spi].astype(float).copy();y=dr['y'][rm[lt],rpi].astype(float).copy();o=gx.get(t,-1);lo=gx.get(lt,-1)
    if o>=0:e[o]=0;y[o]=0
    if lo>=0:y[lo]=0
    xx.append(feature(g,e,ds['profile'],dr['profile'],t));yy.append(y-e);ww.append(1/max(len(train) if split(t)!=0 else len(valid),1));own.append(o);roles.append(split(t)==0);routes.append(s+'->'+r);targets.append(t)
 assert xx
 x=np.asarray(xx);y=np.asarray(yy);fit=~np.array(roles);w=np.asarray(ww);w[fit]/=w[fit].mean();w[~fit]/=w[~fit].mean();mu=x[fit,1:].mean(0);sd=np.maximum(x[fit,1:].std(0),.1);x[:,1:]=(x[:,1:]-mu)/sd
 return dict(x=x,y=y,weight=w,own=np.array(own),validation=~fit,routes=np.array(routes),targets=np.array(targets),feature_mean=mu,feature_scale=sd)

def masked_error(pred,y,own,w):
 mask=np.ones(y.shape,bool);i=np.flatnonzero(own>=0);mask[i,own[i]]=False;e=(pred-y)*mask
 return float(np.sum(w[:,None]*e*e)/np.sum(w[:,None]*mask))

def fit(excluded,shuffle=False):
 g=fit_geometry(excluded);d=paired_table(g,shuffle);tr=~d['validation'];va=~tr;models={};scores={}
 for lam in [1.,10.]:
  coef=old.masked_bilinear(d['x'][tr],g['basis'],d['y'][tr]/g['scale'],d['own'][tr],d['weight'][tr],lam*tr.sum()/d['x'].shape[1]);pred=d['x'][va]@coef@g['basis'].T*g['scale'];scores[str(lam)]=masked_error(pred,d['y'][va],d['own'][va],d['weight'][va]);models[lam]=coef
 lam=min(models,key=lambda v:(scores[str(v)],v));g.update(coef=models[lam],feature_mean=d['feature_mean'],feature_scale=d['feature_scale'],regularization=lam,validation_scores=scores,fit_targets=d['targets'][tr].tolist(),validation_targets=d['targets'][va].tolist(),fit_routes=d['routes'][tr].tolist(),paired_rows=int(tr.sum()),validation_rows=int(va.sum()),shuffled=shuffle)
 assert not set(g['fit_targets'])&set(g['validation_targets']);return g

def estimate(g,query):
 data=old.read_cache();num=np.zeros((300,len(g['axis'])));den=np.zeros_like(num);routes=[];own={str(v):i for i,v in enumerate(data['genes'][g['axis']])}
 from analyze_target_rank_correlations import pairwise_spearman
 for s in g['sources']:
  ds=sources()[s];mp={str(t):i for i,t in enumerate(ds['targets'])};gp={int(v):i for i,v in enumerate(ds['positions'])};p=[gp[int(v)] for v in g['axis']];rho,_,status=pairwise_spearman(query[g['context_axis']],ds['profile'][g['context_axis']],min_pairs=100);weight=np.exp(2*(float(rho) if status=='ok' else -1));family_n=sum(family(t)==family(s) for t in g['sources']);weight/=family_n;features=[];effects=[];ids=[]
  for ti,t in enumerate(data['targets']):
   if str(t) not in mp:continue
   e=ds['y'][mp[str(t)],p].astype(float).copy()
   if str(t) in own:e[own[str(t)]]=0
   features.append(feature(g,e,ds['profile'],query,t));effects.append(e);ids.append(ti)
  if not ids:continue
  x=np.array(features);x[:,1:]=(x[:,1:]-g['feature_mean'])/g['feature_scale'];correction=predict_delta(g,x);vals=np.clip(np.array(effects)+correction,-3,3);num[ids]+=weight*vals;den[ids]+=weight;routes.append({'source':s,'targets':len(ids),'weight':float(weight)})
 return np.divide(num,den,out=np.zeros_like(num),where=den>0),den>0,routes

def predict_delta(g,x):return x@g['coef']@g['basis'].T*g['scale']
def predict(g,query,alpha,parent=None):
 data=old.read_cache();base=old.baseline(data,query,old.excluded_legacy(g['excluded'])) if parent is None else parent.copy();est,mask,routes=estimate(g,query);full=base.copy();full[:,g['axis']]=np.where(mask,est,base[:,g['axis']]);own=np.array([np.flatnonzero(data['genes']==t)[0] for t in data['targets']]);full=np.array([old.gated(v.copy(),np.nan_to_num(query),int(o)) for v,o in zip(full,own)]);pred=np.clip(base+alpha*np.clip(full-base,-1,1),-3,3);pred[np.arange(300),own]=0;return pred,base,routes

def save_model(g,path):
 arrays={k:v for k,v in g.items() if isinstance(v,np.ndarray)};arrays['coef']=g['coef'];np.savez_compressed(path,**arrays);record={k:v for k,v in g.items() if not isinstance(v,np.ndarray) and k!='embedding'};record.update(sha256=sha(path),plan_sha256=sha(A/'plan.json'));write(path.with_suffix('.json'),record)

def freeze():
 assert not (A/'plan.json').exists();deps=['v32_methods.py','run_v32_experiment.py','test_v32_transfer.py','v30_methods.py','run_v30_experiment.py','fit_v7_response.py','v25_methods.py','v24_methods.py','v16_methods.py','v8_methods.py'];inputs={str(p.relative_to(ROOT)):sha(p) for p in (old.A/'sources').glob('*.npz')};inputs['experiments/v7_multisource_adaptive/analysis/source_cache.npz']=sha(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz');inputs['experiments/v24_cd4_response/analysis/models/teacher_all.npz']=sha(ROOT/'experiments/v24_cd4_response/analysis/models/teacher_all.npz');write(A/'plan.json',{'version':32,'dependencies':{p:sha(ROOT/p) for p in deps},'inputs':inputs,'outer_families':OUTERS,'alphas':ALPHAS,'architecture':'Shared16 response programs on<=2048 training-selected shared genes; paired measured intervention tokens; rank-control context coordinates and target STRING features; context-by-response and context-by-target interactions; exact own-label masked bilinear ridge. This adapts the paired-supervision idea, not direct recipient-anchored PerturbMap.','split':'Whole H1, whole Orion HCT116+HEK293T, whole CD4 excluded. Nested alpha selection excludes BOTH outer and inner. Ridge selected on globally identity-held20% paired public training targets; query source response tokens allowed. No fitting to recipient treatment or recipient anchors.','caps':{'basis_rows_per_source':512,'fit_pairs_per_route':256,'validation_pairs_per_route':64},'controls':['paired_identity_shuffled','source_copy','V24'],'gate':{'MSE_ratio':.99,'MAE_ratio':.995,'worst_family_MSE_ratio':1.02,'beats_identity_shuffle_on_both':True},'scope':'Reused public summary LFC development, not blind official count or private anchor scores. No reserve opened; no package/upload. Each later model is tested regardless of this result.','python':sys.version})
def check():
 p=json.loads((A/'plan.json').read_text())
 for n,h in p['dependencies'].items():assert sha(ROOT/n)==h,n
 for n,h in p['inputs'].items():assert sha(ROOT/n)==h,n
 return p