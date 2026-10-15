"""V31: guide-resolved CD4 response teachers, using public raw D2 counts only."""
from pathlib import Path
import json
from functools import lru_cache
import numpy as np
from scipy import sparse,stats
import build_submission as b
from fit_v7_response import read_cache as _read_cache
from v25_methods import raw,baseline
from analyze_target_rank_correlations import pairwise_spearman
ROOT=b.ROOT;A=ROOT/'experiments/v31_guide_response/analysis';CONDITIONS=['CD4_Rest_D2','CD4_Stim48_D2'];ARMS=['guide','pooled','guide_shuffled','target_shuffled']

def write(path,value):
 path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
@lru_cache(maxsize=1)
def data():return _read_cache()
def load_plan():
 p=json.loads((A/'plan.json').read_text())
 for n,h in p['dependencies'].items():assert b.sha256(ROOT/n)==h
 return p

def aggregate(norm,ids):
 x=norm[ids];mu=np.asarray(x.mean(0)).ravel();var=np.maximum(np.asarray(x.power(2).mean(0)).ravel()-mu*mu,0)/max(len(ids)-1,1);return mu,var

def groups(norm,meta,measured,shuffle=False,pool=False):
 labels=meta.gene_target.to_numpy(str);lane=meta['sample'].to_numpy(str);guide=meta.guide_id.to_numpy(str).copy();controls={};control_ids=[]
 for l in sorted(set(lane[labels=='Non-Targeting'])):
  ids=np.flatnonzero((labels=='Non-Targeting')&(lane==l))
  if len(ids)>=20:controls[l]=aggregate(norm,ids);control_ids.extend(ids.tolist())
 assert len(control_ids)>=200
 if shuffle:
  rng=np.random.default_rng(31031)
  for t in sorted(set(labels)-{'Non-Targeting'}):
   for l in sorted(set(lane[labels==t])):
    ids=np.flatnonzero((labels==t)&(lane==l));guide[ids]=guide[rng.permutation(ids)]
 if pool:guide=np.full(len(guide),'pooled',dtype='U6')
 gp={g:i for i,g in enumerate(data()['genes'])};records={};control_profile=np.asarray(norm[control_ids].mean(0)).ravel();control_profile[~measured]=np.nan
 for target in sorted(set(labels)-{'Non-Targeting'}):
  if target not in set(data()['targets']):continue
  block=[]
  for g in sorted(set(guide[labels==target])):
   if g in ['','nan','None','<NA>']:continue
   ids=np.flatnonzero((labels==target)&(guide==g)&np.isin(lane,list(controls)))
   if len(ids)<30:continue
   mu,var=aggregate(norm,ids);base=np.zeros(norm.shape[1]);cv=np.zeros_like(base)
   for l in sorted(set(lane[ids])):
    weight=float((lane[ids]==l).mean());base+=weight*controls[l][0];cv+=weight**2*controls[l][1]
   effect=np.clip(np.log2((mu+.1)/(base+.1)),-3,3);se2=(var/(mu+.1)**2+cv/(base+.1)**2)/np.log(2)**2;se2=np.maximum(se2,.05**2);own=gp[target];kd=float(np.clip(1-(mu[own]+.1)/(base[own]+.1),0,.99)) if base[own]>5 else None;valid=measured&(base>5);valid[own]=False
   block.append({'target':target,'guide':str(g),'n':len(ids),'rows':ids,'effect':effect,'variance':se2,'valid':valid,'knockdown':kd,'baseline':base,'own':own})
  if block:records[target]=block
 return records,control_profile,control_ids

def teacher(blocks,kind):
 effect=np.array([v['effect'] for v in blocks]);variance=np.array([v['variance'] for v in blocks]);valid=np.array([v['valid'] for v in blocks]);n=np.array([v['n'] for v in blocks],float);kd=np.array([.9 if v['knockdown'] is None else v['knockdown'] for v in blocks]);gain=np.clip(.9/np.maximum(kd,.35),.8,1.8) if kind=='guide' else np.ones(len(blocks));effect*=gain[:,None];variance*=gain[:,None]**2
 if kind=='pooled':
  w=n[:,None]*valid;den=w.sum(0);wn=np.divide(w,den,out=np.zeros_like(w),where=den>0);mu=(wn*effect).sum(0);var=(wn*wn*variance).sum(0)+np.sum(wn*(effect-mu)**2,axis=0)
 else:
  first=n[:,None]*valid;den=first.sum(0);wn=np.divide(first,den,out=np.zeros_like(first),where=den>0);start=(wn*effect).sum(0);tau=np.maximum((wn*(effect-start)**2).sum(0)-(wn*variance).sum(0),0);precision=valid/(variance+tau[None,:]+.05**2);lower=np.min(np.where(valid,precision,np.inf),axis=0);lower[~np.isfinite(lower)]=0;precision=np.minimum(precision,10*lower[None,:]);den=precision.sum(0);wn=np.divide(precision,den,out=np.zeros_like(precision),where=den>0);mu=(wn*effect).sum(0);var=np.divide(1,den,out=np.full_like(den,1e6),where=den>0)+tau
 confidence=np.clip(2*stats.norm.cdf(abs(mu)/np.sqrt(np.maximum(var,.0025)))-1,0,1)*n.sum()/(n.sum()+200);confidence[den==0]=0;mu[den==0]=0;return mu.astype(np.float32),confidence.astype(np.float32)

def prepare():
 load_plan();genes=data()['genes'];targets=data()['targets'];out={};calibration=[];metadata=[]
 for ci,condition in enumerate(CONDITIONS):
  x,meta,measured,source=raw(condition,genes);libs=x.sum(1).A1;assert (libs>0).all();norm=x.astype(float).multiply((1e6/libs)[:,None]).tocsr();real,profile,controls=groups(norm,meta,measured);shuffled,_,_=groups(norm,meta,measured,True);pooled,_,_=groups(norm,meta,measured,pool=True);out[f'profile_{ci}']=profile
  for kind in ['guide','pooled','guide_shuffled']:
   mu=np.zeros((300,len(genes)),np.float32);confidence=np.zeros_like(mu)
   for ti,target in enumerate(targets):
    blocks=(shuffled if kind=='guide_shuffled' else pooled if kind=='pooled' else real).get(target,[])
    if not blocks:continue
    mu[ti],confidence[ti]=teacher(blocks,'pooled' if kind=='pooled' else 'guide')
   out[f'{kind}_effect_{ci}']=mu;out[f'{kind}_confidence_{ci}']=confidence
  parent=baseline(data(),profile,['CD4']);own=np.array([np.flatnonzero(genes==t)[0] for t in targets]);per_condition=[]
  for ti,target in enumerate(targets):
   blocks=real.get(target,[])
   if len(blocks)<2:continue
   # Hold out each entire actual guide identity, across every lane.
   for j,actual in enumerate(blocks):
    fitted=[v for k,v in enumerate(blocks) if k!=j];mu,conf=teacher(fitted,'guide');valid=actual['valid'];base=parent[ti];truth=actual['effect'];den=float(((base[valid]-truth[valid])**2).sum());ma=float(abs(base[valid]-truth[valid]).sum())
    if den<=0 or ma<=0:continue
    for alpha in [.125,.25,.5]:
     pred=base+alpha*conf*np.clip(mu-base,-.5,.5);calibration.append({'condition':condition,'target':str(target),'held_guide':actual['guide'],'held_cells':actual['n'],'fitting_guides':[v['guide'] for v in fitted],'alpha':alpha,'mse_ratio':float(((pred[valid]-truth[valid])**2).sum()/den),'mae_ratio':float(abs(pred[valid]-truth[valid]).sum()/ma)})
  metadata.append({'condition':condition,'raw_source':source,'guide_targets':len(real),'multiple_guide_targets':sum(len(v)>=2 for v in real.values()),'guides':[{'target':t,'guide':v['guide'],'n':v['n'],'knockdown':v['knockdown'],'original_rows':v['rows'].tolist()} for t,vs in real.items() for v in vs],'normalization':'Projected competition-axis cell libraries; lane-matched NTC means and sampling variance. Own source suppression supervises dose, never scored as a downstream prediction.'});print('V31 prepared',condition,len(real),sum(len(v)>=2 for v in real.values()),flush=True)
 scores={alpha:np.mean([np.mean([r['mse_ratio'] for r in calibration if r['alpha']==alpha and r['condition']==c]) for c in CONDITIONS]) for alpha in [.125,.25,.5]};alpha=min(scores,key=lambda a:(scores[a],a));file=A/'teacher.npz';np.savez_compressed(file,genes=genes,targets=targets,**out);write(A/'teacher.json',{'sha256':b.sha256(file),'plan_sha256':b.sha256(A/'plan.json'),'metadata':metadata,'calibration':calibration,'chosen_alpha':alpha,'selection_scores':{str(k):float(v) for k,v in scores.items()},'scope':'Guide-held-out calibration within public D2 conditions; correlated culture and controls, not an independent-study or new-context validation.'})

@lru_cache(maxsize=1)
def load_teacher():
 rec=json.loads((A/'teacher.json').read_text());file=A/'teacher.npz';assert b.sha256(file)==rec['sha256'];assert rec['plan_sha256']==b.sha256(A/'plan.json')
 with np.load(file,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files},rec

def predict(query,excluded,arm):
 t,rec=load_teacher();parent=baseline(data(),query,excluded);own=np.array([np.flatnonzero(data()['genes']==v)[0] for v in data()['targets']]);axis=~np.isin(data()['genes'],data()['targets']);rho=np.array([pairwise_spearman(query[axis],t[f'profile_{i}'][axis],min_pairs=100)[0] for i in range(2)]);other=[]
 for ci,family in [(0,'K562'),(4,'H1'),(5,'HCT116'),(6,'HEK293T')]:
  if family not in excluded:other.append(pairwise_spearman(query[axis],data()['profiles'][ci,axis],min_pairs=100)[0])
 best=float(np.nanmax(other)) if np.isfinite(other).any() else -1.;gate=float(np.clip((np.nanmax(rho)-best)/.1,0,1));rho=np.nan_to_num(rho,nan=-1);weight=np.exp((rho-rho.max())/.05);weight/=weight.sum();kind='guide' if arm=='target_shuffled' else arm;mu=np.zeros_like(parent);conf=np.zeros_like(parent);num=np.zeros_like(parent)
 for i in range(2):
  cf=t[f'{kind}_confidence_{i}'];num+=weight[i]*cf*t[f'{kind}_effect_{i}'];conf+=weight[i]*cf
 mu=np.divide(num,conf,out=np.zeros_like(num),where=conf>0)
 if arm=='target_shuffled':
  supported=np.flatnonzero(conf.sum(1)>0);permutation=np.arange(300);permutation[supported]=np.random.default_rng(31032).permutation(supported);mu=mu[permutation];conf=conf[permutation]
 correction=rec['chosen_alpha']*gate*conf*np.clip(mu-parent,-.5,.5);correction[np.arange(300),own]=0;result=np.clip(parent+correction,-3,3).astype(np.float32);return result,parent,{'compatibility_gate':gate,'CD4_rho':rho.tolist(),'other_rho':[float(v) for v in other],'condition_weights':weight.tolist(),'alpha':rec['chosen_alpha'],'changed_targets':int((abs(correction).sum(1)>0).sum())}
