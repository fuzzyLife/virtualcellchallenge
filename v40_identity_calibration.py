"""V40 sign-preserving response calibration trained for intervention identity.

A bounded diagonal map is fitted using measured interventions, not observational
coexpression. No inference-time query treatment or leaderboard score is read.
"""
import argparse,gc,hashlib,json
from pathlib import Path
import numpy as np
import torch
from torch import nn
from torch.nn import functional as F
from threadpoolctl import threadpool_limits
import v38_response_hurdle as parent
from v8_methods import gated
ROOT=parent.ROOT;OUT=ROOT/'experiments/v40_identity_calibration';A=OUT/'analysis'
data=parent.data;b=parent.b;emit=parent.emit
SEED=40040;ARMS=('rank10','rank25');OUTERS=parent.OUTERS
def read(p):return parent.read(p)
def write(p,v):parent.write(p,v)
def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert parent.b.sha256(ROOT/n)==h,n
 for n,h in p['inputs'].items():assert parent.b.sha256(ROOT/n)==h,n
 return p
def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v40_identity_calibration.py','run_v40_identity_calibration.py','run_v38_standard_counts.py','run_v35_v36.py','run_v37_ensemble.py','v38_response_hurdle.py','v15_ensemble_methods.py','v7_response_methods.py','fit_v7_response.py','v16_methods.py','v12_cpu_metrics.py','v12_methods.py','v8_methods.py']
 inputs=['experiments/v7_multisource_adaptive/analysis/source_cache.npz','experiments/v7_multisource_adaptive/analysis/source_provenance.json','experiments/v34_residual_flow/analysis/plan.json','experiments/v38_response_hurdle/analysis/selection.json','experiments/v39_kolf_support/analysis/selection.json']
 write(A/'plan.json',{'dependencies':{s:parent.b.sha256(ROOT/s) for s in deps},'inputs':{s:parent.b.sha256(ROOT/s) for s in inputs},'hypothesis':'Response identity contrastive loss rather than mean-only shrinkage. Two gene-specific parameters, one per response sign. Strictly positive gains bounded0.75..1.25; zero remains zero; own coordinate remains zero. No target-specific new signature invented.','fit':{'epochs':8,'steps_per_epoch':100,'batch_size':32,'learning_rate':.03,'temperature':.15,'MSE_weight':.15,'log_gain_parameter_L2':.05,'seed':SEED,'held_target_hash':'V38 hash0'},'masks':'All300 panel genes excluded from cosine identity loss. Output-own additionally excluded from MSE. Whole outer study excluded from outcomes; recipient study excluded from donor baseline. Biological families sampled uniformly; contexts within families uniformly. Batch-wide valid gene intersection.','controls':'Separately trained within-recipient target-shuffled positive-label model, same recipe. Diagnostic only.','selection':'Compare rank10/rank25 first on all44 public means and pooled-profile discrimination; select larger improvement in equal-study pooled PDS among meanMSE/MAE<=1.01. Then standard V24 count emitter on first4 cases/context, corrected null identical by construction, original fixed V38 six-metric guards. At least two raw metrics improve and beat shuffled diagnostic in summary identity objective. Additional24 raw count cases and secondseed reserved for locked finalist. No private score guarantee; no automatic upload.','known_outcomes':'V38/V39 mean and full standard development outcomes known: better mean error, worse PDS. Transport rejected on severe treatment fidelity loss. Original experiments unchanged.','KOLF_features_used':False,'literature':['https://arxiv.org/html/2606.27752v1','https://arxiv.org/html/2605.19343v1'],'official_score_unknown':True})
def tables(outer,device):
 d=parent.data();panel=np.isin(d['genes'],d['targets']);groups=[]
 for ci in range(7):
  rows=np.flatnonzero(d['row_contexts']==ci)
  if not len(rows):continue
  family=parent.GROUP[str(d['row_families'][rows[0]])]
  if outer!='all' and family==outer:continue
  excluded=tuple(sorted(set(([outer] if outer!='all' else [])+[family])))
  pool=parent.bank(excluded)[-1];own=parent.bank(excluded)[4];query=np.nan_to_num(d['profiles'][ci]);ti=d['row_targets'][rows]
  base=np.asarray([gated(pool[t].copy(),query,int(own[t])) for t in ti],np.float32)
  valid=d['valid'][rows]&np.isfinite(d['profiles'][ci])[None,:]&(query[None,:]>5);valid[np.arange(len(rows)),own[ti]]=False
  val=np.array([parent.split(d['targets'][t])==0 for t in ti]);truth=d['truth'][rows].copy()
  groups.append({'ci':ci,'family':family,'rows':rows,'targets':ti,'validation':val,'base':torch.tensor(base,device=device),'truth':torch.tensor(truth,device=device),'mask':torch.tensor(valid,device=device),'axis':torch.tensor(~panel,device=device),'query':torch.tensor(query,dtype=torch.float32,device=device)})
  print('V40 table',outer,ci,len(rows),flush=True)
 return groups
def loss(theta,g,ids,shuffle=False):
 base=g['base'][ids];truth=g['truth'][ids];mask=g['mask'][ids].all(0);axis=mask&g['axis'];assert int(axis.sum())>=100
 gain=1+.25*torch.tanh(torch.where(base>=0,theta[0],theta[1]));pred=torch.clamp(base*gain,-3,3)
 p=g['query']*.05
 a=torch.log1p(p*torch.exp2(pred))-torch.log1p(p)
 t=torch.log1p(p*torch.exp2(truth))-torch.log1p(p)
 a=F.normalize(a[:,axis],dim=1,eps=1e-8);t=F.normalize(t[:,axis],dim=1,eps=1e-8)
 logits=a@t.T/.15;labels=torch.arange(len(ids),device=base.device)
 if shuffle:labels=torch.roll(labels,1)
 ce=F.cross_entropy(logits,labels);mse=((pred-truth)**2)[:,mask].mean()/torch.clamp(((base-truth)**2)[:,mask].mean(),min=.0025)
 return ce+.15*mse+.05*(theta*theta).mean(),ce,mse
def evaluate(theta,groups,shuffle=False):
 vals=[];families={}
 with torch.no_grad():
  for g in groups:
   ids=np.flatnonzero(g['validation'])
   if len(ids)<3:continue
   _,ce,mse=loss(theta,g,ids,False);f=g['family'];families.setdefault(f,[]).append((float(ce),float(mse)))
 for f,v in families.items():vals.append(np.mean(v,axis=0))
 v=np.mean(vals,axis=0)
 return {'objective':float(v[0]+.15*v[1]),'cross_entropy':float(v[0]),'mse_ratio':float(v[1]),'families':{f:np.mean(v,axis=0).tolist() for f,v in families.items()}}
def train(outer,shuffle=False):
 check();torch.manual_seed(SEED);np.random.seed(SEED);torch.set_num_threads(2);device='cuda' if torch.cuda.is_available() else 'cpu'
 groups=tables(outer,device);theta=nn.Parameter(torch.zeros((2,len(parent.data()['genes'])),device=device));opt=torch.optim.Adam([theta],lr=.03);rng=np.random.default_rng(SEED);families=sorted(set(g['family'] for g in groups));before=evaluate(theta,groups);best=None;bestscore=float('inf');history=[]
 for epoch in range(8):
  losses=[]
  for step in range(100):
   fam=families[int(rng.integers(len(families)))];gg=[g for g in groups if g['family']==fam];g=gg[int(rng.integers(len(gg)))];ids=np.flatnonzero(~g['validation']);ids=rng.choice(ids,min(32,len(ids)),replace=False)
   ll,ce,mse=loss(theta,g,ids,shuffle);opt.zero_grad(set_to_none=True);ll.backward();nn.utils.clip_grad_norm_([theta],1.);opt.step();losses.append(float(ll.detach()))
  val=evaluate(theta,groups);history.append({'epoch':epoch+1,'loss':float(np.mean(losses)),**val})
  if val['objective']<bestscore:bestscore=val['objective'];best=theta.detach().cpu().numpy().copy()
  print('V40',outer,'shuffled' if shuffle else 'identity',epoch+1,val['objective'],val['mse_ratio'],flush=True)
 folder=A/'models';folder.mkdir(exist_ok=True);path=folder/(outer+('_shuffled' if shuffle else '')+'.npz');assert not path.exists()
 np.savez_compressed(path,theta=best,genes=parent.data()['genes'].astype(str))
 write(path.with_suffix('.json'),{'sha256':parent.b.sha256(path),'plan_sha256':parent.b.sha256(A/'plan.json'),'outer':outer,'target_shuffled_training':shuffle,'before':before,'history':history,'best_objective':bestscore,'excluded_recipient_studies':[outer] if outer!='all' else [],'source_recipients':[{'ci':g['ci'],'family':g['family'],'training_targets':[str(parent.data()['targets'][t]) for t in g['targets'][~g['validation']]],'validation_targets':[str(parent.data()['targets'][t]) for t in g['targets'][g['validation']]]} for g in groups]})
 del groups,theta;gc.collect()
 if torch.cuda.is_available():torch.cuda.empty_cache()
def load(outer,shuffle=False):
 p=A/'models'/(outer+('_shuffled' if shuffle else '')+'.npz');r=read(p.with_suffix('.json'));assert parent.b.sha256(p)==r['sha256'] and parent.b.sha256(A/'plan.json')==r['plan_sha256']
 with np.load(p,allow_pickle=False) as z:
  np.testing.assert_array_equal(z['genes'],parent.data()['genes']);return z['theta'].copy()
def predict(theta,base,alpha):
 gain=1+alpha*np.tanh(np.where(base>=0,theta[0],theta[1]));out=np.clip(base*gain,-3,3).astype(np.float32)
 assert not np.any(out*base<0);assert np.all(out[base==0]==0);return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','train']);p.add_argument('--outer',choices=list(OUTERS)+['all']);p.add_argument('--shuffle',action='store_true');a=p.parse_args()
 with threadpool_limits(limits=2):
  if a.stage=='freeze':freeze()
  else:
   for outer in [a.outer] if a.outer else OUTERS:train(outer,a.shuffle)

