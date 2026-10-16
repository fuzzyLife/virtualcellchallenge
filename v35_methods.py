"""V35: gene-level experimental retrieval and a bounded V24 residual.

Independent small model inspired by retrieval-conditioned perturbation work;
not a reproduction of PT-RAG, STATE or their reported results.
"""
import os
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
import copy, gc, hashlib, json
from functools import lru_cache
from pathlib import Path
import numpy as np
import torch
from torch import nn
from scipy import sparse
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits
import v32_methods as old
import v32_runtime as runtime
from v25_methods import raw
from v6_methods import moments

ROOT=old.ROOT
OUT=ROOT/'experiments/v35_response_retrieval'
A=OUT/'analysis'
P36=ROOT/'experiments/v36_population_views/analysis'
SEED=35035
OUTERS=['H1','Orion','CD4']
K=32
EPOCHS=8
STEPS=160
RAW_SOURCES=['H1','HCT116','HEK293T','CD4_Rest_D2','CD4_Stim48_D2']
CI={'H1':4,'HCT116':5,'HEK293T':6,'CD4_Rest':1,'CD4_Stim8hr':2,'CD4_Stim48hr':3}
RAW_NAME={'H1':'H1','HCT116':'HCT116','HEK293T':'HEK293T','CD4_Rest':'CD4_Rest_D2','CD4_Stim48hr':'CD4_Stim48_D2'}
family=old.family
write=old.write
sha=old.sha

def target_split(t):
 return int(hashlib.sha256(('v35-target:'+str(t)).encode()).hexdigest()[:8],16)%5

def setup():
 runtime.install()
 if not hasattr(setup,'done'):
  old.old.embeddings=lru_cache(maxsize=1)(old.old.embeddings)
  setup.done=True
 torch.set_num_threads(2)
 torch.use_deterministic_algorithms(True)
 torch.backends.cuda.matmul.allow_tf32=False

@lru_cache(None)
def vectors():
 k,e=old.old.embeddings()
 return {str(a):np.asarray(b,np.float32) for a,b in zip(k,e)}

def geometry(outer):
 ss=[s for s in old.sources() if outer=='all' or family(s)!=outer]
 common=np.array(sorted(set.intersection(*[set(old.sources()[s]['positions']) for s in ss])))
 data=old.old.read_cache()
 # All common measured outputs: no low-rank response decoder or outcome selection.
 ca=common[~np.isin(data['genes'][common],data['targets'])]
 cf=np.array([rankdata(np.nan_to_num(old.sources()[s]['profile'][ca]))/len(ca) for s in ss])
 mu=cf.mean(0);_,_,vt=np.linalg.svd(cf-mu,full_matrices=False)
 basis=np.zeros((len(ca),4));q=min(4,len(ss)-1);basis[:,:q]=vt[:q].T
 scale=np.maximum(np.std((cf-mu)@basis,0),.01)
 return {'axis':common,'context_axis':ca,'context_mean':mu,'context_basis':basis,'context_scale':scale,'sources':ss,'outer':outer}

def context(g,p):
 v=rankdata(np.nan_to_num(p[g['context_axis']]))/len(g['context_axis'])
 return np.clip((v-g['context_mean'])@g['context_basis']/g['context_scale'],-4,4).astype(np.float32)

def embedding(t,n=64):
 return vectors().get(str(t),np.zeros(64,np.float32))[:n]

@lru_cache(None)
def source_axis(s,axis_tuple):
 d=old.sources()[s];mp={int(g):i for i,g in enumerate(d['positions'])}
 return d['y'][:,[mp[i] for i in axis_tuple]].astype(np.float32)

@lru_cache(None)
def source_embeddings(s):
 return np.array([embedding(n) for n in old.sources()[s]['targets']])

@lru_cache(None)
def candidates(s,t):
 names=old.sources()[s]['targets'];e=source_embeddings(s);te=embedding(t)
 if np.linalg.norm(te)==0:
  ids=np.flatnonzero(names==t);sim=np.zeros(len(names))
 else:
  sim=e@te;ids=np.argsort(-sim,kind='stable')[:4]
  exact=np.flatnonzero(names==t)
  if len(exact) and exact[0] not in ids:ids=np.r_[exact[0],ids[:3]]
 return ids,sim[ids]

def retrieve(g,target,query,excluded,shuffle=False):
 axis=g['axis'];genes=old.old.read_cache()['genes'][axis];q=context(g,query)
 vals=np.zeros((K,len(axis)),np.float32);mask=np.zeros_like(vals,bool);src=np.zeros_like(vals)
 feat=np.zeros((K,79),np.float32);prior=np.full(K,-30.,np.float32);names=[];j=0
 sources=[s for s in g['sources'] if family(s) not in excluded]
 for s in sources:
  ds=old.sources()[s];ids,sim=candidates(s,target);sc=context(g,ds['profile']);bank=source_axis(s,tuple(axis));nt=len(ids)
  # One family vote, then one source vote, rather than rewarding repeated conditions.
  nf=sum(family(n)==family(s) for n in sources)
  perm=np.random.default_rng(SEED).permutation(len(ds['targets'])) if shuffle else np.arange(len(ds['targets']))
  for i,rho in zip(ids,sim):
   if j==K:break
   label=int(perm[i]);ct=str(ds['targets'][i]);donor=str(ds['targets'][label]);same=float(ct==target)
   vals[j]=bank[label];mask[j]=~np.isin(genes,[target,donor]);vals[j,~mask[j]]=0
   src[j]=np.log1p(np.nan_to_num(ds['profile'][axis]))/10
   feat[j]=np.r_[embedding(target,32),embedding(ct,32),sc,q,q-sc,same,float(rho),float(np.log1p(len(ds['targets']))/10)]
   prior[j]=2*float(rho)+same*np.log(4)-np.log(max(nf*nt,1))-.1*float(np.mean((q-sc)**2))
   names.append({'source':s,'candidate':ct,'response_target':donor,'same_target':bool(same),'source_family':family(s)})
   j+=1
 return dict(value=vals,mask=mask,src=src,token=feat,prior=prior,context=q,target=embedding(target,16),names=names)

class RetrievalResidual(nn.Module):
 def __init__(self):
  super().__init__()
  self.attention=nn.Sequential(nn.Linear(79,32),nn.SiLU(),nn.Linear(32,1))
  self.decoder=nn.Sequential(nn.Linear(51,64),nn.SiLU(),nn.Linear(64,32),nn.SiLU(),nn.Linear(32,2))
  nn.init.zeros_(self.attention[-1].weight);nn.init.zeros_(self.attention[-1].bias)
  nn.init.zeros_(self.decoder[-1].weight);nn.init.zeros_(self.decoder[-1].bias)
 def forward(self,b,gene,static=False):
  logits=b['prior']+(0 if static else self.attention(b['token']).squeeze(-1))
  weights=torch.softmax(logits,dim=1)[:,:,None]*b['mask']
  mass=weights.sum(1);w=weights/torch.clamp(mass[:,None,:],min=1e-12)
  value=b['value'];mu=(w*value).sum(1);sd=torch.sqrt(torch.clamp((w*(value-mu[:,None,:])**2).sum(1),min=1e-8))
  pos=(w*(value>0)).sum(1);neg=(w*(value<0)).sum(1)
  same=b['token'][:,:,-3,None]*w;sm=same.sum(1);sv=(same*value).sum(1)/torch.clamp(sm,min=1e-12)
  sp=(w*b['src']).sum(1);base=b['base'];q=b['qlog'];cv=b['control_var'];drop=b['dropout'];known=b['moments_known']
  metrics=torch.stack([base,mu,sd,pos,neg,sv,sm,q,sp,q-sp,base*mu,cv,drop,mu-base,known],dim=-1)
  nn0=base.shape[0];ng=base.shape[1]
  ge=gene[None,:,:].expand(nn0,-1,-1);te=b['target'][:,None,:].expand(-1,ng,-1);cx=b['context'][:,None,:].expand(-1,ng,-1)
  out=self.decoder(torch.cat([metrics,ge,te,cx],dim=-1))
  gate=torch.clamp((torch.expm1(q*10)-1)/5,0,1)*(mass>0)
  delta=.5*torch.tanh(out[:,:,0])*gate
  dispersion=1.5*torch.tanh(out[:,:,1])*gate
  return base+delta,dispersion,delta


def base_for(query,excluded):
 return runtime.baseline(old.old.read_cache(),query,old.old.excluded_legacy(excluded))

def table(g):
 data=old.old.read_cache();axis=g['axis'];rows=[];records=[]
 for s in g['sources']:
  if s not in CI:continue
  ci=CI[s];ex=sorted(set(([g['outer']] if g['outer']!='all' else [])+[family(s)]))
  query=data['profiles'][ci];base=base_for(query,ex)
  idx=np.flatnonzero(data['row_contexts']==ci)
  for ri in idx:
   ti=int(data['row_targets'][ri]);t=str(data['targets'][ti]);tr=retrieve(g,t,query,ex)
   valid=data['valid'][ri,axis].copy()&(data['genes'][axis]!=t)
   if valid.sum()<100 or not tr['mask'].any():continue
   rows.append(dict(**{k:v for k,v in tr.items() if k!='names'},base=base[ti,axis],qlog=np.log1p(np.nan_to_num(query[axis]))/10,truth=data['truth'][ri,axis],valid=valid,control_var=np.zeros(len(axis),np.float32),dropout=np.zeros(len(axis),np.float32),moments_known=np.zeros(len(axis),np.float32)))
   records.append({'source':s,'family':family(s),'target':t,'source_row':int(ri),'target_index':ti,'validation':target_split(t)==0,'excluded':ex,'tokens':tr['names']})
 assert rows
 keys=rows[0].keys();arr={k:np.array([r[k] for r in rows]) for k in keys}
 arr['gene']=np.array([embedding(t,16) for t in data['genes'][axis]])
 return arr,records

def tensor_batch(a,rows,cols,device,views=None):
 out={}
 for k in ['value','mask','src']:
  out[k]=torch.as_tensor(a[k][rows][:,:,cols],device=device,dtype=torch.float32)
 for k in ['token','prior','context','target']:
  out[k]=torch.as_tensor(a[k][rows],device=device,dtype=torch.float32)
 for k in ['base','qlog','control_var','dropout','moments_known','truth','valid']:
  out[k]=torch.as_tensor(a[k][rows][:,cols],device=device,dtype=torch.float32)
 gene=torch.as_tensor(a['gene'][cols],device=device,dtype=torch.float32)
 return out,gene

def balanced_rows(records,ids,n,rng):
 fam=sorted({records[i]['family'] for i in ids});by={f:[i for i in ids if records[i]['family']==f] for f in fam}
 return np.array([rng.choice(by[fam[int(rng.integers(len(fam)))]]) for _ in range(n)])

def loss_mean(pred,truth,base,valid):
 mask=valid;den=torch.clamp(((base-truth)**2*mask).sum(1),min=.02*mask.sum(1))
 sq=((pred-truth)**2*mask).sum(1)/den
 ab=(abs(pred-truth)*mask).sum(1)/torch.clamp((abs(base-truth)*mask).sum(1),min=.1*mask.sum(1))
 return (sq+.2*ab).mean()

def validate(net,a,records,device):
 ids=np.array([i for i,r in enumerate(records) if r['validation']]);results={}
 with torch.no_grad():
  for f in sorted({records[i]['family'] for i in ids}):
   sub=[i for i in ids if records[i]['family']==f];se=be=0.
   for start in range(0,len(sub),8):
    rr=sub[start:start+8]
    for offset in range(0,len(a['gene']),512):
     cols=np.arange(offset,min(offset+512,len(a['gene'])));b,ge=tensor_batch(a,rr,cols,device);p,_,_=net(b,ge);m=b['valid'];se+=float(((p-b['truth'])**2*m).sum());be+=float(((b['base']-b['truth'])**2*m).sum())
   results[f]=se/max(be,1e-12)
 return float(np.mean(list(results.values()))),results

def save(net,g,records,history,best,version,extra=None):
 folder=(A if version==35 else P36)/'models';folder.mkdir(parents=True,exist_ok=True);p=folder/(g['outer']+'.npz');assert not p.exists()
 arrays={'weight__'+k:v.detach().cpu().numpy() for k,v in net.state_dict().items()};arrays.update({k:v for k,v in g.items() if isinstance(v,np.ndarray)})
 np.savez_compressed(p,**arrays)
 fitting=[r for r in records if not r['validation']];validation=[r for r in records if r['validation']]
 write(p.with_suffix('.json'),{'sha256':sha(p),'plan_sha256':sha(A/'plan.json'),'outer':g['outer'],'sources':g['sources'],'selected_epoch':best,'history':history,'fitting_targets':sorted({r['target'] for r in fitting}),'validation_targets':sorted({r['target'] for r in validation}),'training_records':records,'version':version,'extra':extra or {},'outputs':len(g['axis'])})

def load(outer,version=35):
 p=(A if version==35 else P36)/'models'/(outer+'.npz');r=json.loads(p.with_suffix('.json').read_text());assert sha(p)==r['sha256']
 with np.load(p,allow_pickle=False) as z:
  g={k:z[k].copy() for k in z.files if not k.startswith('weight__')};weights={k[8:]:torch.tensor(z[k]) for k in z.files if k.startswith('weight__')}
 g.update(outer=outer,sources=r['sources']);net=RetrievalResidual();net.load_state_dict(weights);net.eval()
 return net,g,r

def train35(outer):
 setup();check();torch.manual_seed(SEED);torch.cuda.manual_seed_all(SEED);rng=np.random.default_rng(SEED)
 g=geometry(outer);a,rec=table(g);device='cuda' if torch.cuda.is_available() else 'cpu';net=RetrievalResidual().to(device);opt=torch.optim.AdamW(net.parameters(),lr=.0005,weight_decay=.001)
 tr=np.array([i for i,r in enumerate(rec) if not r['validation']]);va=[r['target'] for r in rec if r['validation']];assert not {rec[i]['target'] for i in tr}&set(va)
 history=[];score=float('inf');best=None;weights=None
 for epoch in range(1,EPOCHS+1):
  net.train();losses=[]
  for step in range(STEPS):
   ids=balanced_rows(rec,tr,12,rng);cols=rng.choice(len(g['axis']),min(512,len(g['axis'])),replace=False);b,ge=tensor_batch(a,ids,cols,device);p,_,d=net(b,ge);loss=loss_mean(p,b['truth'],b['base'],b['valid'])+.05*(d*d).mean();opt.zero_grad(set_to_none=True);loss.backward();nn.utils.clip_grad_norm_(net.parameters(),1.);opt.step();losses.append(float(loss.detach()))
  net.eval();val,fam=validate(net,a,rec,device);history.append({'epoch':epoch,'loss':float(np.mean(losses)),'validation_MSE_ratio':val,'families':fam})
  if val<score:score=val;best=epoch;weights=copy.deepcopy(net.state_dict())
  print('V35',outer,'epoch',epoch,'validation',val,flush=True)
 net.load_state_dict(weights);save(net,g,rec,history,best,35);del a,net;gc.collect();torch.cuda.empty_cache()


def predict(net,g,query,excluded,targets=None,base=None,shuffle=False,static=False,control=None):
 setup();data=old.old.read_cache();targets=data['targets'] if targets is None else np.asarray(targets);base=base_for(query,excluded) if base is None else base.copy();axis=g['axis'];rows=[];trace=[]
 for i,t in enumerate(targets):
  tr=retrieve(g,str(t),query,excluded,shuffle);trace.append({'target':str(t),'tokens':tr['names']})
  cv=np.zeros(len(axis),np.float32);dr=cv.copy();kn=cv.copy()
  if control is not None:
   pp,vv,dd=control;cv=np.clip(np.log1p(vv[axis]/np.maximum(pp[axis]**2,1e-14))/5,0,2).astype(np.float32);dr=dd[axis];kn=np.ones(len(axis),np.float32)
  rows.append(dict(**{k:v for k,v in tr.items() if k!='names'},base=base[i,axis],qlog=np.log1p(np.nan_to_num(query[axis]))/10,truth=np.zeros(len(axis)),valid=np.ones(len(axis)),control_var=cv,dropout=dr,moments_known=kn))
 arr={k:np.array([r[k] for r in rows]) for k in rows[0]};arr['gene']=np.array([embedding(t,16) for t in data['genes'][axis]])
 result=base.copy();dispersion=np.zeros_like(base);net=net.to('cuda' if torch.cuda.is_available() else 'cpu');device=next(net.parameters()).device
 with torch.no_grad():
  for start in range(0,len(rows),8):
   ids=np.arange(start,min(start+8,len(rows)))
   for offset in range(0,len(axis),512):
    cols=np.arange(offset,min(offset+512,len(axis)));b,ge=tensor_batch(arr,ids,cols,device);p,v,_=net(b,ge,static);result[np.ix_(ids,axis[cols])]=p.cpu().numpy();dispersion[np.ix_(ids,axis[cols])]=v.cpu().numpy()
 for i,t in enumerate(targets):
  own=np.flatnonzero(data['genes']==t)
  if len(own):result[i,own]=0;dispersion[i,own]=0
 return np.clip(result,-3,3).astype(np.float32),dispersion.astype(np.float32),trace

def control_stats(x):
 depth,p,bulk,var=moments(x);drop=1-x.getnnz(0)/x.shape[0]
 return p,var,drop.astype(np.float32)

def check():
 p=json.loads((A/'plan.json').read_text())
 for name,h in p['dependencies'].items():assert sha(ROOT/name)==h,('Frozen dependency changed',name)
 for name,h in p['inputs'].items():assert sha(ROOT/name)==h,('Frozen input changed',name)
 return p

def freeze():
 A.mkdir(parents=True,exist_ok=True);P36.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v35_methods.py','v36_methods.py','run_v35_v36.py','test_v35_v36.py','v32_methods.py','v32_runtime.py','v30_methods.py','v25_methods.py','v24_methods.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py','build_submission.py']
 inputs={str(p.relative_to(ROOT)):sha(p) for p in (old.old.A/'sources').glob('*.npz')}
 for f in ['experiments/v7_multisource_adaptive/analysis/source_cache.npz','experiments/v24_cd4_response/analysis/models/teacher_all.npz','experiments/v34_residual_flow/analysis/plan.json','reference/v30_structured_response/gene_string_embeddings.v0.3.h5']:
  inputs[f]=sha(ROOT/f)
 plan={'version':[35,36],'dependencies':{n:sha(ROOT/n) for n in deps},'inputs':inputs,'outer_families':OUTERS,'seed':SEED,'epochs':EPOCHS,'steps_per_epoch':STEPS,'retrieval':'Up to4 exact/STRING-neighbor targets per allowed source; max32 tokens. Entire recipient family and outer family excluded. Differentiable context-conditioned attention over measured gene-level effects. Remove each response donor own-gene coordinate. Family/source/candidate-normalized prior.','architecture':'Attention79->32->1; per-gene shared decoder51->64->32->2,16 STRING output-gene and16 target coordinates,4 control rank coordinates. All common measured genes, no selected low-rank response decoder. Zero-initial residual head, bounded +/-0.5 log2FC correction; outside modeled/support coordinates exact V24 fallback.','V35':'Family-balanced summary LFC MSE+0.2 MAE with small residual penalty. Fit target hash80%, source-validation target20%; select among8 epochs using source-validation MSE. Whole-family cross-fitted V24 as baseline for every fitting recipient.','V36':'Warm-start its matching V35 model. Eight pairs of disjoint treated populations and matched-lane control populations per training source target. Population LFC consistency supervision weighted by replicate agreement; log-CV2 response head and source NTC variance/dropout features. Summary and population objectives balanced; no observed treated mean/variance at inference.','negative_controls':['V24','retrieved response identities shuffled within source at inference','fixed retrieval attention at inference'],'counts':'V35 mean, V36 mean with unchanged V24 residual10 emitter, V36 learned gene-wise empirical-residual weights0.05..0.25 (Poisson0.65 unchanged), exact libraries and dual mean moments. Same fixed first4 public cases/context plus5 nulls as V34;200-cell D1 descriptive.','mean_gate':{'MSE_ratio_max':.99,'MAE_ratio_max':.995,'worst_family_MSE_ratio_max':1.02,'beats_shuffle_on_both':True},'count_guards':'Family-balanced raw metrics; expression error no more than2% worse; fidelity drop<=0.02; PDS drop<=0.02; null calls no worse than V24. Do not interpret null fidelity as treatment recovery.','release':'User explicitly requests a new package: select best experimental candidate by equal-study mean MSE, subject to count guards when feasible. If no strict improvement passes, label package exploratory/unconfirmed, preserve V24 as champion, do not call it better and do not open reserve. No upload or push.','validation_scope':'Previously spent public families and cases, exploratory; no official private anchors or promised leaderboard improvement. D3 untouched.','literature':[{'url':'https://arxiv.org/abs/2603.07233','idea':'Learn context-dependent retrieval, not fixed semantic averaging; independently adapted with STRING instead of GenePT.'},{'url':'https://arxiv.org/abs/2609.34391','idea':'Independent population views distinguish repeatable effects from sampling variation; moment-summary adaptation, not full author set encoder or reproduction.'}]}
 write(A/'plan.json',plan);write(P36/'plan_reference.json',{'parent_plan_sha256':sha(A/'plan.json')});print('V35/V36 plan frozen',sha(A/'plan.json'),flush=True)