"""Full-axis inference API for the trained V52 joint response model.

This returns mean effects and rank probabilities, never a submission count matrix.
"""
from pathlib import Path
import numpy as np,torch
import v52_joint_response as m
from run_v52_joint_response import load

def query_state(query,outer='all',recipient_study=None):
 query=np.asarray(query,float);assert query.shape==(18533,)
 excluded=tuple(sorted(({outer} if outer!='all' else set())|({recipient_study} if recipient_study else set())))
 desc=m.old.descriptor(query,excluded);base=m.baseline(m.data(),query,m.old.legacy(excluded),production=True)
 return desc,base,m.donor_rank(outer,recipient_study)

def query_features(state,ti,pos):
 desc,base,donors=state;xx,_=m.old.features(desc,ti,pos);zz=[];n=np.zeros(len(pos));total=n.copy();squares=n.copy()
 for k in range(3):
  if k>=len(donors):zz.extend([np.full(len(pos),np.nan) for _ in range(3)]);continue
  study,val,err,mask=donors[k];ok=mask[ti,pos];v=val[ti,pos]
  zz.extend([np.where(ok,v,np.nan),np.where(ok,np.sqrt(err[ti,pos]),np.nan),ok.astype(float)])
  total+=np.where(ok,v,0);squares+=np.where(ok,v*v,0);n+=ok
 prior=np.divide(total,n,out=np.zeros(len(pos)),where=n>0);sd=np.sqrt(np.maximum(np.divide(squares,n,out=np.zeros(len(pos)),where=n>0)-prior**2,0))
 return np.column_stack([xx,base[ti,pos],*zz,prior,sd,n]).astype(np.float32),np.column_stack([base[ti,pos],prior]).astype(np.float32)

def predict(query,outer='all',recipient_study=None,targets=None,device=None):
 m.check();net,ck=load(outer,'joint');device=device or ('cuda' if torch.cuda.is_available() else 'cpu');net.to(device)
 state=query_state(query,outer,recipient_study);d=m.data();targets=d['targets'] if targets is None else np.asarray(targets)
 result_mean=[];result_rank=[];masks=[];pos=np.flatnonzero(np.isfinite(query)&(query>5))
 with torch.no_grad():
  for t in targets:
   ti=int(np.flatnonzero(d['targets']==t)[0]);xx,prior=query_features(state,ti,pos)
   prediction=net(torch.tensor(m.transform(xx,ck['median'],ck['scale']),device=device),torch.tensor(prior,device=device)).cpu().numpy()
   mean=state[1][ti].copy();rank=np.zeros(len(query),np.float32);mask=np.zeros(len(query),bool)
   mean[pos]=prediction[:,0];rank[pos]=prediction[:,1];mask[pos]=True;own=d['genes']==t;mean[own]=0;rank[own]=0;mask[own]=False
   result_mean.append(mean);result_rank.append(rank);masks.append(mask)
 return dict(mean=np.array(result_mean),rank=np.array(result_rank),valid=np.array(masks),genes=d['genes'],targets=targets,
  scope='Learned mean/rank response only; rank effects are not p-values and output arrays are not upload-ready counts.')
