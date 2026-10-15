"""V26 learned fixed-library populations; V24 responses remain frozen."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy import sparse
from v25_methods import raw,baseline
from v22_methods import write,digest,logcpm
from v6_methods import moments
from v16_methods import emit as old_emit
import build_submission as b
ROOT=b.ROOT;OUT=ROOT/'experiments/v26_count_decoder';A=OUT/'analysis'
FAMILIES=['H1','HCT116','HEK293T','CD4'];SOURCES=['H1','HCT116','HEK293T','CD4_Rest_D2','CD4_Stim48_D2']
NAMES=['v24','empirical','independent','conditional','target_shuffled','state_shuffled']
KINDS=['independent','conditional','state_shuffled']
def family_key(name):return 'CD4' if name.startswith('CD4') else name
def load_plan():
 p=json.loads((A/'plan.json').read_text())
 for name,h in p['dependencies'].items():assert b.sha256(ROOT/name)==h,'Frozen V26 code changed: '+name
 for name,h in p['inputs'].items():assert b.sha256(ROOT/name)==h,'Frozen V26 input changed: '+name
 return p
def geometry(key):
 file=ROOT/'experiments/v25_cross_domain/analysis/models'/f'{key}.npz';rec=json.loads(file.with_suffix('.json').read_text());assert b.sha256(file)==rec['sha256'];assert key=='all' or not any(family_key(f)==key for f in rec['training_contexts'])
 with np.load(file,allow_pickle=False) as z:return {k:z[k].copy() for k in ['positions','scale','basis','genes','targets']}
def descriptor(effects,p,geo):
 pos=geo['positions'];cp=1e4*p[pos];der=np.log(2)*cp/(1+cp);return np.clip((effects[...,pos]*der/geo['scale'])@geo['basis'][:,:16],-6,6).astype(np.float32)
def features(count,base,axis):
 libs=np.asarray(count.sum(1)).ravel();v=count[:,axis].toarray().astype(float)*1e4/libs[:,None];return np.clip(np.log1p(v)-np.log1p(1e4*base[:,axis]),-6,6).astype(np.float32)
def state_values(count,meta,ids,state_axis,basis,scale):
 ids=np.asarray(ids);libs=np.asarray(count[ids].sum(1)).ravel();y=count[ids][:,state_axis].toarray().astype(float)*1e4/libs[:,None];np.log1p(y,out=y);batches=meta['sample'].to_numpy(str)[ids]
 for name in set(batches):
  keep=batches==name;pp=np.asarray(count[ids[keep]].multiply((1/libs[keep])[:,None]).mean(0)).ravel();y[keep]-=np.log1p(1e4*pp[state_axis])
 return np.clip(y/scale@basis,-4,4).astype(np.float32)
def fit_probabilities(initial,libs,p,bulk,steps=40):
 # Two-moment exponential tilting with a per-gene damped Newton update.
 q=np.asarray(initial,float).copy();libs=np.asarray(libs,float);active=(p>0)|(bulk>0);q[:,~active]=0;q[:,active]=np.maximum(q[:,active],1e-14);q/=q.sum(1,keepdims=True);weight=libs/libs.mean();w2=weight*weight
 for i in range(steps):
  mean=q.mean(0);pool=weight@q/len(q);second=w2@q/len(q);r=p-mean;t=bulk-pool
  if abs(r).sum()+abs(t).sum()<1e-8:break
  det=np.maximum(mean*second-pool*pool,1e-24);da=(second*r-pool*t)/det;db=(mean*t-pool*r)/det
  # Equal-depth populations need only one independent mean constraint.
  if np.var(weight)<1e-10:da=np.divide(r,mean,out=np.zeros_like(r),where=mean>0);db=np.zeros_like(db)
  step=np.clip(.5*(da[None,:]+weight[:,None]*db[None,:]),-1,1);q*=np.exp(step);q[:,~active]=0;q/=q.sum(1,keepdims=True)
 return q

def integerize(prob,libs,totals,seed):
 # Allocate each row by largest remainder, then preserve its library while
 # repairing exact reference gene totals. Row choice favours existing support.
 rng=np.random.default_rng(seed);libs=np.asarray(libs,np.int64);e=np.asarray(prob,float)*libs[:,None];x=np.floor(e).astype(np.int32)
 for i in range(len(x)):
  n=int(libs[i]-x[i].sum());frac=e[i]-x[i];tie=rng.uniform(0,1e-10,len(frac));take=np.argpartition(frac+tie,-n)[-n:] if n else [];x[i,take]+=1
 delta=np.asarray(totals,np.int64)-x.sum(0,dtype=np.int64);assert delta.sum()==0;neg=list(np.flatnonzero(delta<0));pos=list(np.flatnonzero(delta>0));ii=jj=0;moved=0
 while ii<len(neg) and jj<len(pos):
  a,bg=neg[ii],pos[jj];need=min(int(-delta[a]),int(delta[bg]));rows=np.flatnonzero(x[:,a]>0);priority=prob[rows,bg]/np.maximum(prob[rows,a],1e-14)+rng.uniform(0,1e-12,len(rows));rows=rows[np.argsort(-priority,kind='stable')];available=x[rows,a].astype(np.int64);amount=np.minimum(available,np.maximum(need-np.r_[0,np.cumsum(available)[:-1]],0));x[rows,a]-=amount.astype(np.int32);x[rows,bg]+=amount.astype(np.int32);delta[a]+=need;delta[bg]-=need;moved+=need
  if delta[a]==0:ii+=1
  if delta[bg]==0:jj+=1
 np.testing.assert_array_equal(x.sum(1,dtype=np.int64),libs);np.testing.assert_array_equal(x.sum(0,dtype=np.int64),totals);assert x.min()>=0
 return x,moved

def preserve_depth_margins(sample,reference,seed,bins=16):
 libs=reference.sum(1,dtype=np.int64);order=np.argsort(libs,kind='stable');fixed=np.zeros_like(sample,dtype=np.int32);moved=0
 for i,rows in enumerate(np.array_split(order,min(bins,len(order)))):
  q=np.asarray(sample[rows],float)/libs[rows,None];part,n=integerize(q,libs[rows],reference[rows].sum(0,dtype=np.int64),seed+i);fixed[rows]=part;moved+=n
 np.testing.assert_array_equal(fixed.sum(1,dtype=np.int64),libs);np.testing.assert_array_equal(fixed.sum(0,dtype=np.int64),reference.sum(0,dtype=np.int64))
 return fixed,moved

def repair_cpm(x,reference,max_fraction=.05):
 libs=reference.sum(1,dtype=np.int64);inv=1/libs;target=reference.T@inv;residual=target-x.T@inv;anchors=np.argsort(reference.sum(0))[-8:];moved=0;budget=int(libs.sum()*max_fraction)
 for g in np.r_[np.setdiff1d(np.arange(x.shape[1]),anchors),np.tile(anchors,4)]:
  for attempt in range(8):
   r=residual[g]
   if abs(r)<1e-7 or moved>=budget:break
   sources=np.flatnonzero(x[:,g]>0)
   if not len(sources):break
   order=np.argsort(inv[sources]);sources=np.unique(np.r_[sources[order[:2]],sources[order[-2:]],sources[np.argmax(x[sources,g])]])
   best=None;bestgain=0.
   for a in anchors:
    if a==g:continue
    destinations=np.flatnonzero(x[:,a]>0)
    if not len(destinations):continue
    for donor in sources:
     span=max(abs(inv[destinations]-inv[donor]).max(),1e-12);amounts=np.unique([1,min(int(x[donor,g]),max(1,int(round(abs(r)/span)))),int(x[donor,g])])
     for amount in amounts:
      dst=destinations[(x[destinations,a]>=amount)&(destinations!=donor)]
      if not len(dst):continue
      desired=inv[donor]+r/amount;receiver=int(dst[np.argmin(abs(inv[dst]-desired))]);d=inv[receiver]-inv[donor]
      if d==0:continue
      k=min(int(x[donor,g]),int(x[receiver,a]),max(0,int(round(r/d))),max(0,(budget-moved)//2));shift=k*d;gain=abs(r)-abs(r-shift)
      if k>0 and gain>bestgain:bestgain=gain;best=(donor,receiver,int(a),k,shift)
   if best is None:break
   donor,receiver,a,k,shift=best;x[donor,g]-=k;x[receiver,g]+=k;x[donor,a]+=k;x[receiver,a]-=k;residual[g]-=shift;residual[a]+=shift;moved+=2*k
 np.testing.assert_array_equal(x.sum(1,dtype=np.int64),libs);np.testing.assert_array_equal(x.sum(0,dtype=np.int64),reference.sum(0,dtype=np.int64));assert x.min()>=0
 return x,moved

def decode(model,z,state,condition):
 # CPU NumPy inference avoids GPU scheduling influencing emitted counts.
 w=model['weights'];v=np.concatenate([z,state,condition],axis=1).astype(np.float32)
 for layer in ['d1','d2']:
  v=v@w[layer+'.weight'].T+w[layer+'.bias'];v=v/(1+np.exp(-np.clip(v,-30,30)))
 return (4*np.tanh((v@w['out.weight'].T+w['out.bias'])/4)).astype(float)
def load_decoder(key,kind):
 file=A/'models'/f'{key}_{kind}.npz';rec=json.loads(file.with_suffix('.json').read_text());assert rec['sha256']==b.sha256(file) and rec['plan_sha256']==b.sha256(A/'plan.json')
 with np.load(file,allow_pickle=False) as z:return {'weights':{n[3:]:z[n].copy() for n in z.files if n.startswith('w__')},'output_mask':z['output_mask'].copy(),'state_axis':z['state_axis'].copy(),'state_basis':z['state_basis'].copy(),'state_scale':z['state_scale'].copy(),'geometry':geometry(key),'record':rec}
def learned_initial(template,state,effect,p,model,seed,condition_override=None):
 rng=np.random.default_rng(np.random.SeedSequence([seed,26026]));z=rng.normal(size=(len(state),16)).astype(np.float32);condition=descriptor(effect,p,model['geometry']) if condition_override is None else condition_override;condition=np.repeat(condition[None,:],len(state),axis=0);delta=decode(model,z,state,condition);delta[:,~model['output_mask']]=0
 return delta

def emit(template,effect,own,seed,kind,prepared,atlas,model=None,state=None,condition_override=None,reference=None):
 libs,p,bulk=prepared
 if reference is None:reference,_=old_emit(template,np.arange(len(libs)),effect,own,seed,'residual10',prepared,atlas)
 if kind=='v24':return reference,{'mean_reference_L1':0.,'bulk_exact':True}
 ref=np.asarray(reference.toarray(),np.int32);rp=(ref/libs[:,None]).mean(0);rb=ref.sum(0,dtype=np.int64)/libs.sum();factor=np.exp2(np.clip(effect,-3,3))
 if own is not None:factor[own]=.1
 if kind=='empirical':
  initial=template.toarray().astype(float)*factor;initial/=initial.sum(1,keepdims=True)
 else:
  ss=np.zeros_like(state) if kind=='independent' else state;delta=learned_initial(template,ss,effect,p,model,seed,condition_override);initial=rp[None,:]*np.exp(delta)
  # Coordinates absent in training retain empirical query variation.
  fallback=~model['output_mask'];emp=template[:,fallback].toarray().astype(float)/libs[:,None];initial[:,fallback]=rp[fallback][None,:]*np.divide(emp,p[fallback][None,:],out=np.ones_like(emp),where=p[fallback][None,:]>0)
  if own is not None:initial[:,own]=ref[:,own]/libs
  initial/=np.maximum(initial.sum(1,keepdims=True),1e-15)
 prob=fit_probabilities(initial,libs,rp,rb);rng=np.random.default_rng(seed+260260);sample=np.array([rng.multinomial(int(n),q) for n,q in zip(libs,prob)],dtype=np.int32);fixed,moved=preserve_depth_margins(sample,ref,seed+260260);fixed,cpm_moved=repair_cpm(fixed,ref);mean_error=float(abs((fixed/libs[:,None]).mean(0)-rp).sum());diag={'mean_reference_L1':mean_error,'bulk_exact':True,'depth_strata':16,'bulk_repair_molecules':moved,'CPM_repair_molecules':cpm_moved,'mean_probability_L1':float(abs((fixed/libs[:,None]).mean(0)-p*factor/(p*factor).sum()).sum()),'bulk_probability_L1':float(abs(rb-bulk*factor/(bulk*factor).sum()).sum())}
 # Failure is reported by the stage gate, not repaired with true treatment data.
 return sparse.csr_matrix(fixed),diag
