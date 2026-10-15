"""Masked STRING-bilinear response models and construct-aware source extraction.

Independent SLIM-inspired implementation: source-specific response bases, exact
own-target exclusion in the bilinear loss, and whole-family exclusion at pooling.
"""
from pathlib import Path
import json, hashlib
from functools import lru_cache
import numpy as np
import h5py
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, cg
from sklearn.utils.extmath import randomized_svd
import build_submission as b
from fit_v7_response import read_cache as _read_cache
from v25_methods import baseline
from v8_methods import gated
ROOT=b.ROOT; A=ROOT/'experiments/v30_structured_response/analysis'; R=ROOT/'reference/v30_structured_response'
SOURCES=['K562_full','RPE1_full','H1','HCT116','HEK293T','CD4_Rest','CD4_Stim8hr','CD4_Stim48hr']
FAMILY={'K562_full':'Replogle','RPE1_full':'Replogle','H1':'H1','HCT116':'Orion','HEK293T':'Orion','CD4_Rest':'CD4','CD4_Stim8hr':'CD4','CD4_Stim48hr':'CD4'}
OUTERS=['H1','Orion','CD4']; CONFIGS=[(r,l,a) for r in [16,32] for l in [1.,10.] for a in [.125,.25]]
VARIANTS=['string_full','string_panel','random_full','construct_full']
SEED=30030

@lru_cache(maxsize=1)
def read_cache(): return _read_cache()

def write(path,value):
 path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def load_plan():
 p=json.loads((A/'plan.json').read_text())
 for name,h in p['dependencies'].items():assert b.sha256(ROOT/name)==h,'Frozen code changed: '+name
 return p

def decode(f,node):
 x=f[node]
 if 'categories' in x.attrs:return f[x.attrs['categories']].asstr()[:][x[:]]
 return x.asstr()[:]

def original_source(name,genes):
 rel='reference/K562_gwps_raw_bulk_01.h5ad' if name=='K562_full' else 'reference/v13_resources/rpe1_raw_bulk_01.h5ad'
 with h5py.File(ROOT/rel,'r') as f:
  labels=decode(f,'obs/gene_transcript');syms=decode(f,'var/gene_name');n=f['obs/num_cells_filtered'][:].astype(float);fold=f['obs/fold_expr'][:].astype(float);core=f['obs/core_control'][:].astype(bool);x=f['X'][:].astype(np.float64)
 assert np.isfinite(x).all() and (x>=0).all()
 target=np.array([v.split('_',2)[1] if v.split('_',1)[0].isdigit() and len(v.split('_',2))==3 else '' for v in labels]);controls=core&(target=='non-targeting')&np.isfinite(n)&(n>=30);assert controls.any()
 gp={g:i for i,g in enumerate(genes)};shared=np.array(sorted(set(syms)&set(genes)));positions=np.array([gp[g] for g in shared]);sp={g:i for i,g in enumerate(shared)};sel=np.array([g in sp for g in syms]);mapping=sparse.csr_matrix((np.ones(sel.sum()),(np.flatnonzero(sel),[sp[g] for g in syms[sel]])),shape=(len(syms),len(shared)));x=np.asarray((mapping.T@x.T).T);ctrl=np.average(x[controls],axis=0,weights=n[controls]);control=ctrl*1e6/ctrl.sum();libraries=x.sum(1);base=np.full(len(genes),np.nan);base[positions]=control;usable=(control>5);positions=positions[usable];ctrl=ctrl[usable];control=control[usable];x=x[:,usable]
 valid=np.isfinite(n)&(n>=30)&np.isfinite(fold)&(fold<=.5)&(target!='non-targeting')&(target!='');names=np.array(sorted(set(target[valid])));normal=[];hier=[];total=[];constructs=[]
 for t in names:
  ids=np.flatnonzero(valid&(target==t));xx=x[ids];nn=n[ids];mean=np.average(xx,axis=0,weights=nn);cpm=mean*1e6/np.average(libraries[ids],weights=nn);y=np.clip(np.log2((cpm+.1)/(control+.1)),-3,3);normal.append(y);total.append(nn.sum());constructs.append(len(ids))
  if len(ids)<2:hier.append(y);continue
  cc=xx*1e6/np.maximum(libraries[ids,None],1e-12);effect=np.log2((cc+.1)/(control+.1));gain=np.clip(.9/np.maximum(1-fold[ids],.35),.8,1.8);effect=np.clip(effect*gain[:,None],-3,3)
  # Poisson uncertainty is a bounded proxy; between-construct heterogeneity
  # is retained. Distinct promoter/transcript constructs are not technical replicates.
  se2=(1/np.maximum(xx*nn[:,None],1)+1/np.maximum(ctrl*float(n[controls].sum()),1))/np.log(2)**2*gain[:,None]**2
  w=nn/nn.sum();mu=w@effect;tau=np.maximum(w@((effect-mu)**2)-w@se2,0);prec=1/(se2+tau[None,:]+.05**2);prec=np.minimum(prec,10*prec.min(0,keepdims=True));hier.append((prec*effect).sum(0)/prec.sum(0))
 return {'targets':names,'positions':positions,'y':np.asarray(normal,np.float32),'construct_y':np.asarray(hier,np.float32),'n':np.asarray(total,np.float32),'profile':base,'constructs':np.asarray(constructs,int),'original_labels':labels,'original_target':target,'valid_original_rows':np.flatnonzero(valid),'source_sha256':b.sha256(ROOT/rel)}

def prepare():
 load_plan();data=read_cache();genes=data['genes'];records=[];folder=A/'sources';folder.mkdir(exist_ok=True)
 for name in SOURCES:
  path=folder/(name+'.npz')
  if path.exists():assert b.sha256(path)==json.loads(path.with_suffix('.json').read_text())['sha256'];continue
  if name.endswith('_full'):d=original_source(name,genes);sha=d.pop('source_sha256')
  else:
   ci={'H1':4,'HCT116':5,'HEK293T':6,'CD4_Rest':1,'CD4_Stim8hr':2,'CD4_Stim48hr':3}[name];rows=np.flatnonzero(data['row_contexts']==ci);profile=data['profiles'][ci];pos=np.flatnonzero(np.isfinite(profile)&(profile>5));rows=rows[data['valid'][rows][:,pos].sum(1)>=100];y=data['truth'][rows][:,pos].copy();own=np.array([np.flatnonzero(genes==t)[0] for t in data['targets'][data['row_targets'][rows]]]);mask=data['valid'][rows][:,pos].copy();targetpos={int(v):i for i,v in enumerate(pos)}
   # Require all retained outputs to be observed for each fitting target,
   # allowing only its own-target exclusion. No unmeasured-zero labels.
   for i,o in enumerate(own):
    if int(o) in targetpos:mask[i,targetpos[int(o)]]=True
   pos=pos[mask.all(0)];y=data['truth'][rows][:,pos].copy();prof=profile.copy();d={'targets':data['targets'][data['row_targets'][rows]],'positions':pos,'y':y,'construct_y':y.copy(),'n':np.ones(len(rows),np.float32)*100,'profile':prof,'constructs':np.ones(len(rows),int)};sha=b.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')
  assert len(set(d['targets']))==len(d['targets']) and np.isfinite(d['y']).all();np.savez_compressed(path,**d);write(path.with_suffix('.json'),{'name':name,'family':FAMILY[name],'sha256':b.sha256(path),'source_sha256':sha,'targets':len(d['targets']),'panel_targets':int(np.isin(d['targets'],data['targets']).sum()),'output_genes':len(d['positions']),'multiple_effective_construct_targets':int((d['constructs']>=2).sum()),'scope':'Pseudobulk means are not integer single cells. Missing outputs absent. Whole-CD4 and whole-Orion exclusions enforced at inference. Own target excluded in fitted objective.'});print('V30 source',name,len(d['targets']),len(d['positions']),flush=True)
 records=[json.loads((folder/(n+'.json')).read_text()) for n in SOURCES];write(A/'source_inventory.json',records)

def embeddings():
 with h5py.File(R/'gene_string_embeddings.v0.3.h5','r') as f:
  keys=np.array(sorted(f.keys()));vectors=np.array([f[k][:64] for k in keys],float)
 assert vectors.shape[1]==64 and np.isfinite(vectors).all();vectors/=np.maximum(np.linalg.norm(vectors,axis=1,keepdims=True),1e-12)
 return keys,vectors

def masked_bilinear(x,g,y,own,weight,penalty):
 """Solve weighted bilinear ridge with each row's own-gene label absent."""
 x=np.asarray(x,float);g=np.asarray(g,float);y=np.asarray(y,float);own=np.asarray(own,int);w=np.asarray(weight,float);gram=x.T@(w[:,None]*x);gg=g.T@g;rhs=x.T@(w[:,None]*(y@g));ok=own>=0;xx=x[ok];go=g[own[ok]];ww=w[ok]
 if len(xx):rhs-=xx.T@((ww*y[np.flatnonzero(ok),own[ok]])[:,None]*go)
 # y[own] is zero-filled only as storage; correction removes that coordinate's
 # model term from the normal matrix, so zero is NEVER a supervised own label.
 def mv(flat):
  m=flat.reshape(x.shape[1],g.shape[1]);z=gram@m@gg+penalty*m
  if len(xx):z-=xx.T@((ww*np.einsum('ij,ij->i',xx@m,go))[:,None]*go)
  return z.ravel()
 op=LinearOperator((rhs.size,rhs.size),matvec=mv,dtype=np.float64);sol,info=cg(op,rhs.ravel(),rtol=1e-8,atol=1e-10,maxiter=400);assert info==0,('CG failed',info);assert np.linalg.norm(mv(sol)-rhs.ravel())/max(np.linalg.norm(rhs),1e-12)<1e-6
 return sol.reshape(rhs.shape)

def fit_one(name,variant):
 load_plan();data=read_cache();keys,e=embeddings();lookup={k:i for i,k in enumerate(keys)};d=np.load(A/'sources'/f'{name}.npz',allow_pickle=False);ts=d['targets'];keep=np.array([t in lookup for t in ts]);keep&=np.isin(ts,data['targets']) if variant=='string_panel' else True;ids=np.flatnonzero(keep)
 if len(ids)<20:
  file=A/'models'/f'{name}_{variant}.npz';file.parent.mkdir(exist_ok=True);empty={f'r{r}_l{int(l)}':np.zeros((300,len(d['positions'])),np.float32) for r in [16,32] for l in [1.,10.]};np.savez_compressed(file,positions=d['positions'],profile=d['profile'],present=np.zeros(300,bool),targets=data['targets'],**empty);write(file.with_suffix('.json'),{'sha256':b.sha256(file),'source_sha256':b.sha256(A/'sources'/f'{name}.npz'),'plan_sha256':b.sha256(A/'plan.json'),'family':FAMILY[name],'source':name,'variant':variant,'training_targets':ts[ids].tolist(),'training_target_count':len(ids),'auxiliary_targets':0,'abstained':True,'reason':'Fewer than20 eligible source targets with embeddings'});print('V30 abstained',name,variant,len(ids),flush=True);return
 ts=ts[ids];pos=d['positions'];y=d['construct_y' if variant=='construct_full' else 'y'][ids].astype(float);w=d['n'][ids].astype(float);w=w/(w+100);w/=w.mean();gx={str(g):i for i,g in enumerate(data['genes'][pos])};own=np.array([gx.get(t,-1) for t in ts]);mask=np.ones(y.shape,bool);ii=np.flatnonzero(own>=0);mask[ii,own[ii]]=False;y[~mask]=0
 bias=(w[:,None]*y).sum(0)/(w[:,None]*mask).sum(0);centered=y-bias;centered[~mask]=0;scale=np.maximum(np.sqrt((w[:,None]*centered**2).sum(0)/(w[:,None]*mask).sum(0)),.1);scaled=centered/scale
 rank=min(32,len(ts)-2,len(pos)-2);_,_,vt=randomized_svd(scaled*np.sqrt(w[:,None]),n_components=rank,random_state=SEED,n_iter=4);basis=vt.T
 ee=e.copy()
 if variant=='random_full':ee=ee[np.random.default_rng(SEED).permutation(len(ee))]
 x=ee[[lookup[t] for t in ts]];mu=x.mean(0);sd=np.maximum(x.std(0),.01);x=(x-mu)/sd;query=data['targets'];present=np.array([t in lookup for t in query]);qx=np.zeros((len(query),64));qx[present]=(ee[[lookup[t] for t in query[present]]]-mu)/sd
 outputs={};models={}
 for rr in [16,32]:
  r=min(rr,rank);g=basis[:,:r]
  for lam in [1.,10.]:
   ww=masked_bilinear(x,g,scaled,own,w,lam*len(ts)/64);pred=((qx@ww)@g.T)*scale+bias;pred=np.clip(pred,-3,3);pred[~present]=0
   for ti,t in enumerate(query):
    if t in gx:pred[ti,gx[t]]=0
   key=f'r{rr}_l{int(lam)}';outputs[key]=pred.astype(np.float32);models[key]=ww
 file=A/'models'/f'{name}_{variant}.npz';file.parent.mkdir(exist_ok=True);np.savez_compressed(file,positions=pos,profile=d['profile'],present=present,targets=query,basis=basis,scale=scale,bias=bias,feature_mean=mu,feature_scale=sd,**{'weight_'+k:v for k,v in models.items()},**outputs);write(file.with_suffix('.json'),{'sha256':b.sha256(file),'source_sha256':b.sha256(A/'sources'/f'{name}.npz'),'plan_sha256':b.sha256(A/'plan.json'),'family':FAMILY[name],'source':name,'variant':variant,'training_targets':ts.tolist(),'training_target_count':len(ts),'auxiliary_targets':int((~np.isin(ts,query)).sum()),'own_missing_labels':int((own>=0).sum()),'basis_scope':'Training source only. Own coordinates filled with centered zero for SVD initialization; exact own-coordinate exclusion in ridge normal equations. No held-family source used at prediction.'});print('V30 fitted',name,variant,len(ts),len(pos),flush=True)

@lru_cache(maxsize=64)
def model_arrays(name,variant):
 file=A/'models'/f'{name}_{variant}.npz';rec=json.loads(file.with_suffix('.json').read_text());assert b.sha256(file)==rec['sha256']
 with np.load(file,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files if not k.startswith('weight_') and k not in ['basis','scale','bias','feature_mean','feature_scale']}

def predict_graph(query,excluded,variant,rank,lam):
 data=read_cache();num=np.zeros((300,len(data['genes'])));den=np.zeros_like(num);details=[];models=[]
 from analyze_target_rank_correlations import pairwise_spearman
 panel=~np.isin(data['genes'],data['targets'])
 for name in SOURCES:
  if FAMILY[name] in excluded:continue
  file=A/'models'/f'{name}_{variant}.npz';rec=json.loads(file.with_suffix('.json').read_text());assert rec['family'] not in excluded and b.sha256(file)==rec['sha256'];z=model_arrays(name,variant);rho,_,status=pairwise_spearman(query[panel],z['profile'][panel],min_pairs=100);rho=float(rho) if status=='ok' else -1.;models.append((name,z,rho))
 if not models:return num.astype(np.float32),den>0,details
 best=max(v[2] for v in models)
 families={f:np.zeros(num.shape,np.uint8) for f in set(FAMILY[n] for n,_,_ in models)}
 for name,z,rho in models:
  p=z['positions'];valid=(np.nan_to_num(query[p])>1)[None,:]&z['present'][:,None];families[FAMILY[name]][:,p]+=valid.astype(np.uint8)
 for name,z,rho in models:
  p=z['positions'];q=np.nan_to_num(query[p]);valid=(q>1)[None,:]&z['present'][:,None];weight=np.exp(2*(rho-best));coord_weight=weight/np.maximum(families[FAMILY[name]][:,p],1);value=z[f'r{rank}_l{int(lam)}'];num[:,p]+=value*valid*coord_weight;den[:,p]+=valid*coord_weight;details.append({'source':name,'family':FAMILY[name],'raw_weight':float(weight),'rho':rho,'normalization':'Divide by available source count within family separately per target/gene'})
 return np.divide(num,den,out=np.zeros_like(num),where=den>0).astype(np.float32),den>0,details

def excluded_legacy(excluded):
 result=[]
 for key in excluded:result+= {'Orion':['HCT116','HEK293T'],'Replogle':['K562'],'H1':['H1'],'CD4':['CD4']}[key]
 return result

def predict(query,excluded,variant,config,parent=None):
 data=read_cache();base=baseline(data,query,excluded_legacy(excluded)) if parent is None else parent.copy();r,lam,alpha=config;est,mask,details=predict_graph(query,excluded,variant,r,lam);own=np.array([np.flatnonzero(data['genes']==t)[0] for t in data['targets']]);est=np.array([gated(v.copy(),np.nan_to_num(query),int(o)) for v,o in zip(est,own)])
 correction=np.clip(est-base,-.5,.5)*mask;correction[np.arange(300),own]=0;result=np.clip(base+alpha*correction,-3,3).astype(np.float32);result[np.arange(300),own]=0
 return result,base,{'sources':details,'changed_targets':int((abs(result-base).sum(1)>0).sum()),'config':[r,lam,alpha]}
