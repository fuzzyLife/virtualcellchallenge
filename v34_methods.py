"""V34 trainable-mean control-anchored residual flow in measured RNA space.

An independently implemented conditional-flow pilot. No observed treated mean,
query response descriptor, or treated posterior is supplied at inference.
"""
import os
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8')
from pathlib import Path
import json,copy,gc,hashlib,time
import numpy as np
from scipy import sparse
from scipy.optimize import linear_sum_assignment
import torch
from torch import nn
from threadpoolctl import threadpool_limits
import v32_methods as m
from v25_methods import raw
ROOT=m.ROOT;A=ROOT/'experiments/v34_residual_flow/analysis';SEED=34034
SOURCES=['H1','HCT116','HEK293T','CD4_Rest_D2','CD4_Stim48_D2']
def family(s):return 'Orion' if s in ['HCT116','HEK293T'] else 'CD4' if s.startswith('CD4') else s

def target_feature(t):
 if not hasattr(target_feature,'table'):
  keys,e=m.old.embeddings();target_feature.table={str(k):v[:32].astype(np.float32) for k,v in zip(keys,e)}
 return target_feature.table.get(str(t),np.zeros(32,np.float32))
from functools import lru_cache
target_feature=lru_cache(None)(target_feature)

def logspace(x,axis):
 libs=x.sum(1).A1;assert (libs>0).all();return np.log1p(x[:,axis].toarray().astype(np.float64)*1e4/libs[:,None]).astype(np.float32)
def normmean(x):return np.asarray(x.astype(float).multiply((1e6/x.sum(1).A1)[:,None]).mean(0)).ravel()

class Flow(nn.Module):
 def __init__(self,genes,conditions):
  super().__init__();self.body=nn.Sequential(nn.Linear(2*genes+conditions+1,128),nn.SiLU(),nn.Linear(128,128),nn.SiLU(),nn.Linear(128,genes));nn.init.zeros_(self.body[-1].weight);nn.init.zeros_(self.body[-1].bias)
 def forward(self,x,x0,t,c):return self.body(torch.cat([x,x0,t,c],dim=1))

def integrate(net,x0,c,steps=8):
 x=x0.clone()
 for i in range(steps):
  t=torch.full((len(x),1),i/steps,device=x.device);x=x+net(x,x0,t,c)/steps
 return x

def coupled(x,y,own,return_order=False):
 xx=x.copy();yy=y.copy()
 if own>=0:xx[:,own]=0;yy[:,own]=0
 cost=np.maximum(np.sum(xx*xx,1)[:,None]+np.sum(yy*yy,1)[None,:]-2*xx@yy.T,0);_,order=linear_sum_assignment(cost);y=y[order].copy()
 if own>=0:y[:,own]=x[:,own]
 return (y,order) if return_order else y

def prepare(outer):
 g=m.fit_geometry([outer]);axis=g['axis'][-1024:];d=m.old.read_cache();genes=d['genes'];gp={str(t):i for i,t in enumerate(genes[axis])};records=[];groups=[];source_controls=[];centers=[];cells=[];offset=0;context_profiles=[];source_names=[];originals=[]
 for si,s in enumerate([s for s in SOURCES if family(s)!=outer]):
  x,meta,measured,rec=raw(s,genes);assert measured[axis].all(),('Unmeasured model output',s);labels=meta.gene_target.to_numpy(str);lanes=meta['sample'].to_numpy(str);rng=np.random.default_rng(SEED+si);control_ids=[];lane_center={};lane_pool={}
  for lane in sorted(set(lanes[labels=='Non-Targeting'])):
   ids=np.flatnonzero((labels=='Non-Targeting')&(lanes==lane));ids=rng.permutation(ids)[:64]
   if len(ids)<8:continue
   lane_pool[lane]=ids;lane_center[lane]=np.log1p(normmean(x[ids])[axis]*.01).astype(np.float32);control_ids.extend(ids.tolist())
  assert control_ids;profile=normmean(x[control_ids]);profile[~measured]=np.nan;context_profiles.append(profile);context=m.ctx(g,profile).astype(np.float32);source_names.append(s);source_controls.append(np.concatenate([logspace(x[ids],axis)-lane_center[l] for l,ids in lane_pool.items()]));ctl=source_controls[-1];selected=[];source_groups=[];source_centers=[]
  for target in sorted(set(labels)-{'Non-Targeting'}):
   if target not in set(d['targets']):continue
   ids=np.flatnonzero((labels==target)&np.isin(lanes,list(lane_pool)));ids=rng.permutation(ids)[:48]
   if len(ids)<16:continue
   z=logspace(x[ids],axis)-np.array([lane_center[l] for l in lanes[ids]]);condition=np.r_[target_feature(target),context].astype(np.float32);group={'source':s,'family':family(s),'source_index':si,'target':str(target),'validation':m.split(target)==0,'null':False,'own':gp.get(str(target),-1),'condition':condition,'y':z,'centers':np.array([lane_center[l] for l in lanes[ids]]),'original_rows':ids.tolist()};groups.append(group);selected.extend(ids.tolist());source_groups.append(str(target))
  # Null groups use empirical source controls; zero target descriptor.
  groups.append({'source':s,'family':family(s),'source_index':si,'target':'Non-Targeting','validation':False,'null':True,'own':-1,'condition':np.r_[np.zeros(32),context].astype(np.float32),'y':ctl,'centers':np.array([np.mean(list(lane_center.values()),0)]*len(ctl)),'original_rows':control_ids})
  records.append({'source':s,'family':family(s),'source_record':rec,'control_rows':control_ids,'treatment_rows':selected,'targets':source_groups});originals.append(rec);del x;gc.collect();print('V34 prepared',outer,s,'groups',len(source_groups),'controls',len(ctl),flush=True)
 folder=A/'training';folder.mkdir(exist_ok=True);file=folder/f'{outer}.npz';arrays={'axis':axis,'context_axis':g['context_axis'],'context_mean':g['context_mean'],'context_basis':g['context_basis'],'context_scale':g['context_scale'],'source_profiles':np.array(context_profiles),'source_names':np.array(source_names)};group_records=[]
 for i,v in enumerate(groups):
  arrays[f'y_{i}']=v['y'];arrays[f'center_{i}']=v['centers'];arrays[f'condition_{i}']=v['condition'];group_records.append({k:v for k,v in v.items() if k not in ['y','centers','condition']})
 for i,v in enumerate(source_controls):arrays[f'control_{i}']=v
 np.savez_compressed(file,**arrays);record={'sha256':m.sha(file),'plan_sha256':m.sha(A/'plan.json'),'outer':outer,'source_records':records,'groups':group_records,'axis_genes':genes[axis].tolist(),'basis_training':g['basis_training'],'normalization':'Per-cell full aligned-axis library, log1p(CP10k), subtract matched-lane NTC log mean. Own target suppressed label removed from coupling and model inputs/loss.'};m.write(file.with_suffix('.json'),record);return arrays,record

def load_training(outer):
 file=A/'training'/f'{outer}.npz'
 if not file.exists():return prepare(outer)
 rec=json.loads(file.with_suffix('.json').read_text());assert m.sha(file)==rec['sha256'] and rec['plan_sha256']==m.sha(A/'plan.json')
 with np.load(file,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files},rec

def train(outer):
 plan=check();a,rec=load_training(outer);groups=rec['groups'];fit=[i for i,v in enumerate(groups) if not v['validation']];valid=[i for i,v in enumerate(groups) if v['validation']];assert valid and not ({groups[i]['target'] for i in fit}&{groups[i]['target'] for i in valid});torch.manual_seed(SEED);torch.cuda.manual_seed_all(SEED);torch.set_num_threads(2);torch.use_deterministic_algorithms(True);torch.backends.cuda.matmul.allow_tf32=False;device='cuda' if torch.cuda.is_available() else 'cpu';net=Flow(len(a['axis']),len(a['condition_0'])).to(device);opt=torch.optim.AdamW(net.parameters(),lr=.001,weight_decay=.0001);rng=np.random.default_rng(SEED);families=sorted({groups[i]['family'] for i in fit});byfam={f:[i for i in fit if groups[i]['family']==f] for f in families};history=[];best=float('inf');bestweights=None
 def batch(i,validation=False):
  v=groups[i];rr=np.random.default_rng(SEED+i) if validation else rng;n=min(32,len(a[f'y_{i}']));yi=rr.choice(len(a[f'y_{i}']),n,replace=False);ctl=a[f'control_{v["source_index"]}'];x0=ctl[rr.choice(len(ctl),n,replace=len(ctl)<n)].copy();y,order=coupled(x0,a[f'y_{i}'][yi],v['own'],True);yi=yi[order];c=np.repeat(a[f'condition_{i}'][None,:],n,0);mask=np.ones(len(a['axis']),np.float32)
  if v['own']>=0:mask[v['own']]=0
  return tuple(torch.tensor(z,dtype=torch.float32,device=device) for z in [x0,y,c,mask,a[f'center_{i}'][yi]])
 for epoch in range(8):
  net.train();losses=[];start=time.time()
  for step in range(400):
   f=families[step%len(families)];i=int(rng.choice(byfam[f]));x0,y,c,mask,center=batch(i);t=torch.tensor(rng.uniform(0,1,(len(x0),1)),dtype=torch.float32,device=device);xt=x0+(y-x0)*t;velocity=net(xt,x0,t,c);truth=y-x0;flow=((velocity-truth)**2*mask).sum()/mask.sum()/len(x0);vm=(velocity*mask).mean(0);tm=(truth*mask).mean(0);mean=((vm-tm)**2).mean();direction=1-torch.nn.functional.cosine_similarity(vm[None,:],tm[None,:]).mean();loss=flow+mean+(0 if groups[i]['null'] else .1)*direction;opt.zero_grad(set_to_none=True);loss.backward();torch.nn.utils.clip_grad_norm_(net.parameters(),5.);opt.step();assert torch.isfinite(loss);losses.append(float(loss.detach()))
  net.eval();scores=[]
  with torch.no_grad():
   for i in valid:
    x0,y,c,mask,center=batch(i,True);pred=integrate(net,x0,c);anchor=torch.tensor(np.log1p(np.nan_to_num(a['source_profiles'][groups[i]['source_index']][a['axis']])*.01),dtype=torch.float32,device=device);pm=torch.expm1(torch.clamp(pred+anchor,0,12)).mean(0);ym=torch.expm1(torch.clamp(y+center,0,12)).mean(0);base=torch.expm1(torch.clamp(x0+anchor,0,12)).mean(0);err=((pm-ym)**2*mask).sum();den=((base-ym)**2*mask).sum();scores.append((groups[i]['family'],float(err/max(float(den),1e-8))))
  value=float(np.mean([np.mean([v for f,v in scores if f==ff]) for ff in sorted(set(f for f,v in scores))]));history.append({'epoch':epoch+1,'training_loss':float(np.mean(losses)),'held_target_inference_mean_error_ratio_vs_controls':value,'seconds':time.time()-start});print('V34 epoch',outer,epoch+1,history[-1],flush=True)
  if value<best:best=value;bestweights={k:v.detach().cpu().numpy().copy() for k,v in net.state_dict().items()}
 folder=A/'models';folder.mkdir(exist_ok=True);file=folder/f'{outer}.npz';np.savez_compressed(file,**{'w__'+k:v for k,v in bestweights.items()},**{k:a[k] for k in ['axis','context_axis','context_mean','context_basis','context_scale']});m.write(file.with_suffix('.json'),{'sha256':m.sha(file),'plan_sha256':m.sha(A/'plan.json'),'training_sha256':rec['sha256'],'outer':outer,'fitting_sources':[v['source'] for v in rec['source_records']],'fitting_targets':sorted({groups[i]['target'] for i in fit}),'validation_targets':sorted({groups[i]['target'] for i in valid}),'history':history,'selected_epoch':min(history,key=lambda r:r['held_target_inference_mean_error_ratio_vs_controls'])['epoch'],'selection':'Eight fixed epochs, minimum actual8-step generated held-target source mean-CP10k error. No outer-treatment outcomes or true means supplied to predictor.','GPU':torch.cuda.get_device_name(0),'torch':torch.__version__});del net,a;gc.collect();torch.cuda.empty_cache()

def load_model(outer):
 file=A/'models'/f'{outer}.npz';rec=json.loads(file.with_suffix('.json').read_text());assert m.sha(file)==rec['sha256']
 with np.load(file,allow_pickle=False) as z:a={k:z[k].copy() for k in z.files}
 net=Flow(len(a['axis']),32+a['context_basis'].shape[1]);net.load_state_dict({k[3:]:torch.tensor(v) for k,v in a.items() if k.startswith('w__')});net.eval();return net,a,rec

def generated(net,a,template,query,target=None,shuffle=False):
 axis=a['axis'];center=np.log1p(np.nan_to_num(query[axis])*.01).astype(np.float32);x0=logspace(template,axis)-center;context=m.ctx(a,query).astype(np.float32);emb=np.zeros(32,np.float32) if target is None else target_feature(str(target));
 if shuffle and target is not None:
  ts=m.old.read_cache()['targets'];perm=np.random.default_rng(SEED).permutation(len(ts));ti=int(np.flatnonzero(ts==target)[0]);emb=target_feature(str(ts[perm[ti]]))
 c=np.repeat(np.r_[emb,context][None,:],len(x0),0).astype(np.float32);pieces=[]
 with torch.no_grad():
  for lo in range(0,len(x0),128):pieces.append(integrate(net,torch.tensor(x0[lo:lo+128]),torch.tensor(c[lo:lo+128])).numpy())
 values=np.expm1(np.clip(np.concatenate(pieces)+center,0,12))*100.;assert np.isfinite(values).all();return values

def check():
 p=json.loads((A/'plan.json').read_text())
 for n,h in p['dependencies'].items():assert m.sha(ROOT/n)==h,n
 return p

def freeze():
 assert not (A/'plan.json').exists();old=json.loads((ROOT/'experiments/v25_cross_domain/analysis/plan.json').read_text());contexts=old['contexts']
 for i,c in enumerate(contexts):
  c['excluded']=['HCT116','HEK293T'] if family(c['family'])=='Orion' else [family(c['family'])];rr=np.random.default_rng(SEED+i).permutation(c['left']+c['right']);c['null']={'target':'Non-Targeting','seed':SEED+i,'actual':rr[:400].tolist(),'reference':rr[400:800].tolist(),'template':rr[800:1200].tolist(),'fit_controls':rr[800:].tolist()}
 deps=['v34_methods.py','run_v34_experiment.py','test_v34_flow.py','v32_methods.py','v32_runtime.py','v25_methods.py','v12_cpu_metrics.py','v12_methods.py','v8_methods.py','v16_methods.py','v6_methods.py'];m.write(A/'plan.json',{'version':34,'dependencies':{n:m.sha(ROOT/n) for n in deps},'source_inputs':old['source_inputs'],'contexts':contexts,'sources':SOURCES,'architecture':{'gene_outputs':1024,'hidden':[128,128],'condition':'32 STRING target coordinates +4-or-fewer training-control rank PCs','epochs':8,'steps_per_epoch':400,'batch':32,'integration_steps':8,'optimizer':'AdamW lr0.001 wd0.0001','loss':'Own-masked OT-coupled velocity MSE +1 mean-velocity MSE +0.1 direction cosine (non-null only)'},'split':'Entire H1, entire Orion and entire CD4 excluded. Target-identity20% source validation using V32 hash; no source fitting response for those identities. Source controls reusable. Query conditioning only target identity/embedding and NTC context, consistent with fitting; no observed mean or response descriptor.','mean_arms':['v24','flow_mean25','flow_mean100','target_shuffle25'],'primary':'flow_mean25','count_arms':['v24','flow_mean25','flow_shape25'],'mean_gate':'Primary MSE<=0.99 and MAE<=0.995, worst-family MSE<=1.02, CPM MSE<=1.02; beat target shuffle on both LFC errors. Count inference is tested regardless on a fixed diagnostic subset; failure cannot advance to release.','counts':'All frozen public mean cases; first4 cases/context count diagnostic with actual400 for H1/Orion/Stim and200 for D1. One400-cell null/context. Identical baseline emitter for learned mean; shape25 blends neural profiles with baseline counts before exact-library largest-remainder rounding. Raw Wilcoxon/BH and panel metrics, no private anchors.','raw_guards':old['raw_guards'],'limitations':'Partial1024-gene pilot; outside modeled axis V24 fallback. K562/RPE1 have pseudobulk inputs here, not raw flow training cells. Matched-lane centering with source controls; OT coupling across centered populations is a modeling device, not observed individual counterfactuals. Reused public biological cases; no new blind lineage.','release':'No reserve or package/upload during these tests. Test this model independently of V32/V33 outcomes.'})