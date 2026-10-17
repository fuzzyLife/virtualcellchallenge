"""V52 joint mean/rank response learning; new files only, D3 never accessed."""
import os
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
from pathlib import Path
from functools import lru_cache
import gc,hashlib,json,platform,importlib.metadata
import numpy as np
from scipy import sparse,stats
import torch
from torch import nn
import build_submission as b
import v38_response_hurdle as old
from v25_methods import raw,baseline
ROOT=b.ROOT;A=ROOT/'experiments/v52_joint_mean_rank/analysis'
SEED=52052
SOURCES=('H1','HCT116','HEK293T','CD4_Rest_D2','CD4_Stim48_D2')
GROUP={'H1':'H1','HCT116':'Orion','HEK293T':'Orion','CD4_Rest_D2':'CD4','CD4_Stim48_D2':'CD4'}
OUTERS=('H1','Orion','CD4')
ARMS=('joint','mean_only','rank_only','shuffled_joint')
CONFIG=dict(hidden=[128,64],epochs=8,steps_per_epoch=192,batch_size=2048,learning_rate=.0005,
 weight_decay=.001,dropout=.05,pairs_per_group=768,treated_per_lane_view=64,
 controls_per_lane_view=128,minimum_treated_per_view=16,minimum_controls_per_lane_view=16,
 mean_residual_bound=2.,rank_residual_bound=.25)

def write(p,value):b.write_json(p,value)
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def target_split(t):return int(hashlib.sha256(('v52-target:'+str(t)).encode()).hexdigest()[:8],16)%5

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v52_joint_response.py','run_v52_joint_response.py','test_v52_joint_response.py','v38_response_hurdle.py','v25_methods.py','v24_methods.py','v16_methods.py','v15_ensemble_methods.py','fit_v7_response.py']
 inputs=['experiments/v7_multisource_adaptive/analysis/source_provenance.json','experiments/v24_cd4_response/analysis/models/teacher_all.json']
 write(A/'plan.json',dict(version=52,config=CONFIG,sources=SOURCES,groups=GROUP,arms=ARMS,seed=SEED,
  dependencies={s:b.sha256(ROOT/s) for s in deps},inputs={s:b.sha256(ROOT/s) for s in inputs},
  scope='Shared mean/rank predictor on raw source populations; full study held out. Previously examined public biology, not blind validation. No decoder or submission claimed.',
  labels='Two disjoint treated/control views within each lane. CPM mean log2 ratio with .1 pseudocount, clipped +/-3; theta=U/(n_t*n_c)-.5 with half credit for ties. Lane weights follow treated counts. Output gene must have both control-view means>5 CPM and be measured. Own gene excluded.',
  features='Recipient controls and old donor features; exclude whole outer and recipient studies from every donor. Raw donor rank effects are pooled with one vote per biological study. Missing values have explicit indicators. No recipient response labels in features.',
  validation='All hash0 targets excluded from network/ridge fitting, retained for source inner validation. Allowed other-study responses to those targets remain retrieval features: cross-context, not unseen-intervention evaluation. Hold out HCT116+HEK293T together and every CD4 source together.',
  selection='Fixed eight epochs; select epoch only by source hash0 joint normalized mean/rank MSE (single-head uses its head loss). Report each outer separately, equal-study summaries, linear baseline, trained shuffled-label control and view reliability. Never select a release from these metrics alone.',
  shuffle='Within each source, cyclic target-label permutation restricted to fitting targets. Same output gene; inputs unchanged. Validation and outer test always genuine labels.',
  pilot_success='Mean MSE no worse than V24; rank MSE at least 5% below pooled rank, below zero rank, and both MSEs below trained shuffle. Not sufficient for submission release.',
  D3_opened=False,automatic_upload=False,python=platform.python_version(),versions={n:importlib.metadata.version(n) for n in ['numpy','scipy','torch','scikit-learn']}))
 print('V52 frozen',b.sha256(A/'plan.json'),flush=True)

def check():
 p=read(A/'plan.json');assert p['config']==CONFIG
 for n,h in p['dependencies'].items():assert b.sha256(ROOT/n)==h,n
 for n,h in p['inputs'].items():assert b.sha256(ROOT/n)==h,n
 return p

@lru_cache(1)
def data():return old.data()

def cpm(x,rows):
 z=x[np.asarray(rows)].astype(float);lib=z.sum(1).A1;assert (lib>0).all()
 return z.multiply((1e6/lib)[:,None]).tocsr()

def rank_effect(x,c):
 """Columns are genes; ties receive half credit, without expression pseudocounts."""
 x=np.asarray(x,float);c=np.asarray(c,float);n,m=len(x),len(c)
 ranks=stats.rankdata(np.concatenate([x,c]),axis=0,method='average')
 return (ranks[:n].sum(0)-n*(n+1)/2)/(n*m)-.5

def summarize_views(x,views,measured,own):
 means=[];ranks=[];controls=[];totals=[]
 for lanes in views:
  total=sum(len(t) for t,c in lanes);mu=np.zeros(x.shape[1]);cm=mu.copy();theta=mu.copy()
  for tr,cr in lanes:
   tx=cpm(x,tr);cx=cpm(x,cr);weight=len(tr)/total
   tm=tx.mean(0).A1;cc=cx.mean(0).A1;mu+=weight*tm;cm+=weight*cc
   cols=np.flatnonzero(measured)
   # Compute rank effects on all measured genes that can enter the view's gate.
   for st in range(0,len(cols),256):
    pos=cols[st:st+256];theta[pos]+=weight*rank_effect(tx[:,pos].toarray(),cx[:,pos].toarray())
  controls.append(cm);means.append(np.clip(np.log2((mu+.1)/(cm+.1)),-3,3));ranks.append(theta);totals.append(total)
 valid=measured&(np.min(controls,axis=0)>5);valid[own]=False
 return np.array(means,np.float32),np.array(ranks,np.float32),valid,np.mean(controls,axis=0).astype(np.float32),totals

def prepare(source):
 check();folder=A/'labels';folder.mkdir(exist_ok=True);path=folder/(source+'.npz')
 if path.exists():assert b.sha256(path)==read(path.with_suffix('.json'))['sha256'];return
 d=data();genes=d['genes'];targets=d['targets'];x,meta,measured,receipt=raw(source,genes)
 labels=meta.gene_target.to_numpy(str);lanes=meta['sample'].to_numpy(str)
 pools={s:np.flatnonzero((labels=='Non-Targeting')&(lanes==s)) for s in sorted(set(lanes[labels=='Non-Targeting']))}
 pools={k:v for k,v in pools.items() if len(v)>=2*CONFIG['minimum_controls_per_lane_view']}
 ntc=np.concatenate(list(pools.values()));query=cpm(x,ntc).mean(0).A1;query[~measured]=np.nan
 yy=np.zeros((300,2,len(genes)),np.float32);rr=np.zeros_like(yy);valid=np.zeros((300,len(genes)),bool);queries=np.zeros((300,len(genes)),np.float32);records=[]
 for ti,t in enumerate(targets):
  rng=np.random.default_rng(SEED+ti);views=[[],[]];trace=[]
  for lane,pool in pools.items():
   tr=np.flatnonzero((labels==t)&(lanes==lane));n=min(len(tr)//2,CONFIG['treated_per_lane_view']);nc=min(len(pool)//2,CONFIG['controls_per_lane_view'])
   if n<2:continue
   tr=rng.permutation(tr)[:2*n];cr=rng.permutation(pool)[:2*nc]
   for v in range(2):views[v].append((tr[v*n:(v+1)*n],cr[v*nc:(v+1)*nc]))
   trace.append(dict(lane=lane,treated_views=[tr[:n].tolist(),tr[n:].tolist()],control_views=[cr[:nc].tolist(),cr[nc:].tolist()]))
  if not views[0] or sum(len(t) for t,c in views[0])<CONFIG['minimum_treated_per_view']:continue
  own=int(np.flatnonzero(genes==t)[0]);ym,yr,mask,q,nn=summarize_views(x,views,measured,own)
  yy[ti]=ym;rr[ti]=yr;valid[ti]=mask;queries[ti]=q
  records.append(dict(target=str(t),target_index=ti,cells_per_view=nn,valid_genes=int(mask.sum()),lanes=trace))
  if len(records)%20==0:print('V52 rank labels',source,len(records),'target',t,flush=True)
 np.savez_compressed(path,mean=yy,rank=rr,valid=valid,query=query,target_control=queries,genes=genes,targets=targets)
 write(path.with_suffix('.json'),dict(sha256=b.sha256(path),plan_sha256=b.sha256(A/'plan.json'),source=source,study=GROUP[source],source_record=receipt,groups=records,
  mean_label='Clipped batch-matched CPM log2FC (.1 only for mean labels)',rank_label='Batch-matched U probability minus .5; no rank pseudocount',
  rows=int(x.shape[0]),measured_genes=int(measured.sum()),D3_opened=False))
 print('V52 labels completed',source,len(records),flush=True);del x;gc.collect()

@lru_cache(5)
def labels(source):
 p=A/'labels'/(source+'.npz');r=read(p.with_suffix('.json'));assert r['plan_sha256']==b.sha256(A/'plan.json') and b.sha256(p)==r['sha256']
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def donor_rank(outer,recipient):
 excluded={recipient} if outer=='all' else {outer,recipient};families=[]
 for study in OUTERS:
  if study in excluded:continue
  values=[];masks=[];noise=[]
  for source in SOURCES:
   if GROUP[source]!=study:continue
   z=labels(source);values.append(z['rank'].mean(1));masks.append(z['valid']);noise.append((z['rank'][:,0]-z['rank'][:,1])**2)
  den=np.sum(masks,axis=0);val=np.divide(np.sum(np.where(masks,values,0),axis=0),den,out=np.zeros_like(values[0]),where=den>0)
  err=np.divide(np.sum(np.where(masks,noise,0),axis=0),den,out=np.zeros_like(values[0]),where=den>0)
  families.append((study,val,err,den>0))
 return families

@lru_cache(16)
def context_features(source,outer):
 z=labels(source);query=z['query'];excluded=tuple(sorted({GROUP[source]}|({outer} if outer!='all' else set())))
 desc=old.descriptor(query,excluded)
 legacy=old.legacy(excluded);base=baseline(data(),query,legacy,production=True)
 return desc,base,donor_rank(outer,GROUP[source])

def features(source,outer,ti,pos):
 desc,base,donors=context_features(source,outer);xx,_=old.features(desc,ti,pos)
 zz=[];n=np.zeros(len(pos));total=n.copy();squares=n.copy()
 for k in range(3):
  if k>=len(donors):zz.extend([np.full(len(pos),np.nan) for _ in range(3)]);continue
  family,val,err,mask=donors[k];ok=mask[ti,pos];v=val[ti,pos]
  zz.extend([np.where(ok,v,np.nan),np.where(ok,np.sqrt(err[ti,pos]),np.nan),ok.astype(float)])
  total+=np.where(ok,v,0);squares+=np.where(ok,v*v,0);n+=ok
 prior=np.divide(total,n,out=np.zeros(len(pos)),where=n>0)
 sd=np.sqrt(np.maximum(np.divide(squares,n,out=np.zeros(len(pos)),where=n>0)-prior**2,0))
 xx=np.column_stack([xx,base[ti,pos],*zz,prior,sd,n]).astype(np.float32)
 return xx,np.column_stack([base[ti,pos],prior]).astype(np.float32)

def training_table(outer):
 xs=[];ys=[];ws=[];priors=[];records=[];sources=[];tis=[];positions=[]
 for source in SOURCES:
  if outer!='all' and GROUP[source]==outer:continue
  z=labels(source)
  for ti in np.flatnonzero(z['valid'].any(1)):
   rng=np.random.default_rng(SEED+int(ti));pos=np.flatnonzero(z['valid'][ti]);pos=rng.choice(pos,min(CONFIG['pairs_per_group'],len(pos)),replace=False)
   xx,prior=features(source,outer,int(ti),pos);y=np.column_stack([z['mean'][ti].mean(0)[pos],z['rank'][ti].mean(0)[pos]])
   noise=np.column_stack([(z['mean'][ti,0,pos]-z['mean'][ti,1,pos])**2,(z['rank'][ti,0,pos]-z['rank'][ti,1,pos])**2])
   rel=np.clip(np.array([.1,.01])/(np.array([.1,.01])+noise),.1,1)/len(pos)
   xs.append(xx);ys.append(y);ws.append(rel);priors.append(prior);sources.extend([source]*len(pos));tis.extend([int(ti)]*len(pos));positions.extend(pos.tolist())
   records.append(dict(source=source,study=GROUP[source],target=str(z['targets'][ti]),target_index=int(ti),pairs=len(pos),validation=target_split(z['targets'][ti])==0))
  print('V52 table',outer,source,sum(len(x) for x in xs),flush=True)
 x=np.concatenate(xs);y=np.concatenate(ys).astype(np.float32);w=np.concatenate(ws).astype(np.float32);prior=np.concatenate(priors)
 sources=np.array(sources);ti=np.array(tis);pos=np.array(positions);val=np.array([target_split(data()['targets'][i])==0 for i in ti])
 # Equal source study, then equal source context, then equal target, then equal pair.
 for study in sorted({GROUP[s] for s in sources}):
  ss=sorted({s for s in sources if GROUP[s]==study})
  for s in ss:
   ids=sources==s;nt=len(np.unique(ti[ids]));w[ids]/=len(ss)*nt
 w/=w.mean(0)
 return dict(x=x,y=y,w=w,prior=prior,source=sources,ti=ti,pos=pos,val=val,records=records)

class Network(nn.Module):
 def __init__(self,n):
  super().__init__();self.trunk=nn.Sequential(nn.Linear(n,128),nn.SiLU(),nn.LayerNorm(128),nn.Dropout(.05),nn.Linear(128,64),nn.SiLU());self.mean=nn.Linear(64,1);self.rank=nn.Linear(64,1)
  nn.init.zeros_(self.mean.weight);nn.init.zeros_(self.mean.bias);nn.init.zeros_(self.rank.weight);nn.init.zeros_(self.rank.bias)
 def forward(self,x,prior,arm='joint'):
  h=self.trunk(x);m=2*torch.tanh(self.mean(h).squeeze(1)/2);r=.25*torch.tanh(self.rank(h).squeeze(1))
  if arm=='rank_only':m=m*0
  if arm=='mean_only':r=r*0
  return torch.stack([torch.clamp(prior[:,0]+m,-3,3),torch.clamp(prior[:,1]+r,-.5,.5)],1)

def normalizer(x):
 med=np.nanmedian(x,axis=0);med=np.nan_to_num(med);v=np.where(np.isfinite(x),x,med);scale=np.maximum(np.std(v,axis=0),.05)
 return med.astype(np.float32),scale.astype(np.float32)
def transform(x,med,scale):
 missing=~np.isfinite(x);v=np.where(missing,med,x)
 return np.column_stack([np.clip((v-med)/scale,-8,8),missing.astype(np.float32)]).astype(np.float32)

def shuffled_labels(table):
 out=table['y'].copy();perms={}
 for s in sorted(set(table['source'])):
  ts=np.unique(table['ti'][(table['source']==s)&~table['val']]);perms[s]={int(t):int(v) for t,v in zip(ts,np.roll(ts,1))}
  z=labels(s)
  for ti,dest in perms[s].items():
   ids=np.flatnonzero((table['source']==s)&(table['ti']==ti)&~table['val']);pos=table['pos'][ids]
   ok=z['valid'][dest,pos];out[ids[ok],0]=z['mean'][dest].mean(0)[pos[ok]];out[ids[ok],1]=z['rank'][dest].mean(0)[pos[ok]]
   # Missing shuffled target/gene pairs abstain from the training loss.
   table['shuffle_available'][ids]=ok
 return out,perms
