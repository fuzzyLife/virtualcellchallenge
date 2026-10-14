"""V22 bounded cell-state/continuous-response pilot; no prior code mutated."""
from pathlib import Path
import json,hashlib,gc
import numpy as np,pandas as pd
from scipy import sparse
from sklearn.utils.extmath import randomized_svd
from v6_methods import moments
from v16_methods import emit as v16_emit,source_arrays
from v8_methods import gated
from fit_v7_response import raw_inputs,read_cache
import build_submission as b
ROOT=b.ROOT;OUT=ROOT/'experiments/v22_state_response';A=OUT/'analysis'
FAMILIES=['H1','HCT116','HEK293T']
NAMES=['v16','scalar','independent','conditional_half','conditional','shuffled_state','amplitude']

def clean(x):
 if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [clean(v) for v in x]
 if isinstance(x,np.ndarray):return clean(x.tolist())
 if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
 if isinstance(x,(int,np.integer)):return int(x)
 if isinstance(x,(bool,np.bool_)):return bool(x)
 return x

def write(path,value):b.write_json(path,clean(value))
def digest(x):
 h=hashlib.sha256()
 for v in (x.data,x.indices,x.indptr):h.update(v.tobytes())
 return h.hexdigest()
def logcpm(x):
 libs=np.asarray(x.sum(1)).ravel();assert (libs>0).all()
 y=x.multiply((1e6/libs)[:,None]).tocsr();y.data=np.log1p(y.data);return y

def load_plan():
 p=json.loads((A/'plan.json').read_text())
 for name,h in p['dependencies'].items():assert b.sha256(ROOT/name)==h, 'Frozen V22 dependency changed: '+name
 return p

def raw(family,genes):
 if family!='CD4':return raw_inputs(family,genes)
 p=A/'cd4_selected';rec=json.loads((p/'extraction.json').read_text());assert b.sha256(p/'counts.npz')==rec['counts_sha256']
 x=sparse.load_npz(p/'counts.npz');meta=pd.read_csv(p/'cells.csv');axis=pd.read_csv(p/'genes.csv')
 assert b.sha256(p/'cells.csv')==rec['cells_sha256'];np.testing.assert_array_equal(axis.gene_name.to_numpy(str),genes)
 return x,meta,axis.measured.to_numpy(bool),rec

def fit_model(excluded):
 plan=load_plan();folder=A/'models';folder.mkdir(exist_ok=True);dest=folder/(excluded+'.npz')
 if dest.exists():return load_model(excluded)
 data=read_cache();genes=data['genes'];targets=data['targets'];training=[f for f in FAMILIES if f!=excluded]
 assert excluded not in training
 # Stable state features are selected only from training intervention summaries.
 rows=np.isin(data['row_families'],training);valid=data['valid'][rows];value=np.where(valid,abs(data['truth'][rows]),np.nan)
 with np.errstate(invalid='ignore'):sensitivity=np.nanmedian(value,axis=0)
 sensitivity=np.nan_to_num(sensitivity,nan=10.);nonpanel=~np.isin(genes,targets)
 common=np.isfinite(data['profiles'][np.isin(np.array(['K562','CD4_Rest','CD4_Stim8hr','CD4_Stim48hr','H1','HCT116','HEK293T']),training)]).all(0)
 eligible=common&nonpanel&(np.nanmin(data['profiles'][np.isin(np.array(['K562','CD4_Rest','CD4_Stim8hr','CD4_Stim48hr','H1','HCT116','HEK293T']),training)],axis=0)>5)
 state_axis=np.flatnonzero(eligible);state_axis=state_axis[np.argsort(sensitivity[state_axis],kind='stable')[:768]]
 assert len(state_axis)>=128
 output=np.flatnonzero(common&nonpanel);output=output[np.argsort(np.nanmax(np.where(valid,abs(data['truth'][rows]),0),axis=0)[output],kind='stable')[-2048:]]
 controls=[];records={}
 for fi,f in enumerate(training):
  x,meta,measured,rec=raw(f,genes);y=logcpm(x);labels=meta.gene_target.to_numpy(str);batches=meta['sample'].to_numpy(str);cs=np.flatnonzero(labels=='Non-Targeting');bn=sorted(set(batches));means={};good=[]
  for batch in bn:
   ids=cs[batches[cs]==batch]
   if len(ids)>=8:means[batch]=np.asarray(y[ids][:,state_axis].mean(0)).ravel();good.extend(ids.tolist())
  take=np.random.default_rng(2200+fi).choice(good,min(700,len(good)),replace=False)
  c=y[take][:,state_axis].toarray()-np.array([means[batches[i]] for i in take]);controls.append(c)
  records[f]={'raw':rec,'control_rows':take.tolist(),'good_batches':list(means)}
  del x,y;gc.collect()
 c=np.concatenate(controls);scale=np.maximum(np.std(c,axis=0),.2);_,_,vt=randomized_svd(c/scale,n_components=4,random_state=22022);basis=vt.T
 accum={k:[] for k in ('interaction','strength_beta','strength_residual','strength_cv','support')};diagnostics=[]
 for fi,f in enumerate(training):
  x,meta,measured,rec=raw(f,genes);y=logcpm(x);labels=meta.gene_target.to_numpy(str);batches=meta['sample'].to_numpy(str);cs=np.flatnonzero(labels=='Non-Targeting');batch_parameters={}
  for batch in records[f]['good_batches']:
   ids=cs[batches[cs]==batch];mu=np.asarray(y[ids].mean(0)).ravel();zz=(y[ids][:,state_axis].toarray()-mu[state_axis])/scale@basis;zs=np.maximum(zz.std(0),.2);zz/=zs
   w=np.linalg.solve(zz.T@zz+40*np.eye(4),zz.T@(y[ids][:,output].toarray()-mu[output]))
   batch_parameters[batch]=(mu,zs,w)
  control_var=np.maximum(np.asarray(y[cs][:,output].power(2).mean(0)).ravel()-np.asarray(y[cs][:,output].mean(0)).ravel()**2,.1)
  inter=np.zeros((len(targets),4,len(output)),np.float32);sb=np.zeros((len(targets),4),np.float32);sr=np.zeros((len(targets),128),np.float32);cv=np.zeros(len(targets),np.float32);support=np.zeros(len(targets),np.int32)
  # First collect residuals; the across-target mean is the shared treatment component.
  pieces={}
  for ti,target in enumerate(targets):
   ids=np.flatnonzero((labels==target)&np.isin(batches,list(batch_parameters)))
   if len(ids)<80:continue
   if len(ids)>800:ids=np.random.default_rng(22000+ti).choice(ids,800,replace=False)
   z=np.zeros((len(ids),4));res=np.zeros((len(ids),len(output)))
   for batch in sorted(set(batches[ids])):
    where=np.flatnonzero(batches[ids]==batch);mu,zs,w=batch_parameters[batch];q=y[ids[where]]
    zz=(q[:,state_axis].toarray()-mu[state_axis])/scale@basis/zs;zz=np.clip(zz,-3,3)
    z[where]=zz;res[where]=q[:,output].toarray()-mu[output]-zz@w
   pieces[ti]=(ids,z,res)
  generic=np.mean([r[2].mean(0) for r in pieces.values()],axis=0) if pieces else np.zeros(len(output))
  for ti,(ids,z,res) in pieces.items():
   fingerprint=res.mean(0)-generic;signal=abs(fingerprint)/np.sqrt(control_var)
   axis=np.argsort(signal)[-256:];direction=fingerprint[axis]/control_var[axis];norm=float(direction@fingerprint[axis])
   if norm<=.1:continue
   # Cross-gene agreement shrinks a noisy inferred response strength toward one.
   aa=axis[::2];bb=axis[1::2]
   def projection(j):
    w=fingerprint[j]/control_var[j];return (res[:,j]-generic[j])@w/max(float(w@fingerprint[j]),1e-8)
   sa=projection(aa);sd=projection(bb);cov=max(float(np.cov(sa,sd,ddof=1)[0,1]),0);noise=max(float(np.var(sa-sd)/4),1e-4);reliability=cov/(cov+noise)
   strength=1+reliability*(.5*(sa+sd)-1);strength=np.clip(strength,.15,2.);strength/=strength.mean()
   beta=np.linalg.solve(z.T@z+80*np.eye(4),z.T@(strength-1));residual=strength-1-z@beta
   design=strength[:,None]*z;design-=design.mean(0)
   remainder=res-generic-strength[:,None]*fingerprint
   fitted=np.linalg.solve(design.T@design+160*np.eye(4),design.T@remainder)
   # Output slope converts log1p-CPM residuals to approximate log2-fold changes.
   mu=np.mean([v[0][output] for v in batch_parameters.values()],axis=0);derivative=np.maximum(1-np.exp(-mu),.2)*np.log(2)
   fitted/=derivative[None,:];fitted=np.clip(fitted,-.4,.4)*len(ids)/(len(ids)+200)
   inter[ti]=fitted;sb[ti]=np.clip(beta,-.2,.2);sr[ti]=np.quantile(residual,np.linspace(0,1,128));cv[ti]=strength.var();support[ti]=len(ids)
   diagnostics.append({'family':f,'target':str(targets[ti]),'cells':len(ids),'strength_cv2':float(cv[ti]),'cross_gene_reliability':reliability,'interaction_norm':float(np.linalg.norm(fitted)),'shared_program_norm':float(np.linalg.norm(generic)),'specific_fingerprint_norm':float(np.linalg.norm(fingerprint))})
  for key,val in zip(accum,(inter,sb,sr,cv,support)):accum[key].append(val)
  del x,y,pieces;gc.collect();print('V22 conditional source fit',excluded,f,'supported',int((support>0).sum()),flush=True)
 sup=np.stack(accum['support']);weight=sup/(sup+200);den=weight.sum(0)
 result={'genes':genes,'targets':targets,'state_axis':state_axis,'output_axis':output,'scale':scale,'basis':basis,'support':sup.sum(0)}
 for key in ('interaction','strength_beta','strength_residual','strength_cv'):
  arr=np.stack(accum[key]);w=weight.reshape(weight.shape+(1,)*(arr.ndim-2));d=den.reshape(den.shape+(1,)*(arr.ndim-2));result[key]=np.divide((arr*w).sum(0),d,out=np.zeros(arr.shape[1:]),where=d>0).astype(np.float32)
 np.savez_compressed(dest,**result);write(dest.with_suffix('.json'),{'excluded_family':excluded,'training_families':training,'original_inputs':records,'plan_sha256':b.sha256(A/'plan.json'),'sha256':b.sha256(dest),'state_features':len(state_axis),'output_features':len(output),'supported_targets':int((result['support']>0).sum()),'scope':'Bounded factor pilot; baseline-state proxy can retain treatment effects. No individual causal pairing claim.'});pd.DataFrame(diagnostics).to_csv(folder/(excluded+'_fits.csv'),index=False)
 return result

def load_model(excluded):
 p=A/'models'/(excluded+'.npz');r=json.loads(p.with_suffix('.json').read_text());assert r['sha256']==b.sha256(p) and r['plan_sha256']==b.sha256(A/'plan.json');assert excluded not in r['training_families']
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def query_states(x,model):
 y=logcpm(x)[:,model['state_axis']].toarray();mu=y.mean(0);z=(y-mu)/model['scale']@model['basis'];sd=np.maximum(z.std(0),.2);return np.clip(z/sd,-3,3)

def conditional_delta(z,effect,model,ti,kind,seed,query):
 n=len(z);rng=np.random.default_rng(np.random.SeedSequence([seed,22022]));res=rng.choice(model['strength_residual'][ti],n);zz=z.copy()
 if kind=='shuffled_state':zz=zz[rng.permutation(n)]
 if kind in ('scalar','independent'):strength=1+res
 else:strength=1+res+zz@model['strength_beta'][ti]
 strength=np.clip(strength,.15,2.);strength/=strength.mean();delta=(strength[:,None]-1)*effect[None,:]
 if kind not in ('scalar',):
  if kind=='independent':zz=zz[rng.permutation(n)]
  design=strength[:,None]*zz;design-=design.mean(0);part=design@model['interaction'][ti];gate=np.clip(query[model['output_axis']]/5,0,1);delta[:,model['output_axis']]+=part*gate
 delta-=delta.mean(0);delta=np.clip(delta,-.5,.5)*(.125 if kind=='conditional_half' else .25)
 return delta

def emit(template,effect,own,seed,kind,prepared,atlas,model=None,ti=None,z=None):
 n=template.shape[0]
 if kind=='v16' or model is None or model['support'][ti]==0 or not np.any(effect):return v16_emit(template,np.arange(n),effect,own,seed,'residual10',prepared,atlas)
 if kind=='amplitude':return v16_emit(template,np.arange(n),effect*1.05,own,seed,'residual10',prepared,atlas)
 libs,p,bulk=prepared;factor=np.exp2(np.clip(effect,-3,3))
 if own is not None:factor[own]=.1
 wp=p*factor;wp/=wp.sum();wb=bulk*factor;wb/=wb.sum();delta=conditional_delta(z,effect,model,ti,kind,seed,p*1e6)
 if own is not None:delta[:,own]=0
 cell_factor=factor[None,:]*np.exp2(delta);prob=p[None,:]*cell_factor;prob/=prob.sum(1,keepdims=True);expected=libs[:,None]*prob;empirical=template.toarray().astype(float)*cell_factor;empirical*=np.divide(libs,empirical.sum(1))[:,None]
 initial=.25*expected+.65*np.random.default_rng(seed).poisson(expected)+.10*empirical;initial/=initial.sum(1,keepdims=True)
 count=atlas.dual_moment_counts(initial,wp,wb,depths=np.asarray(libs,np.int64),seed=int(seed),iterations=200)
 np.testing.assert_array_equal(count.sum(1,dtype=np.int64),libs);assert np.issubdtype(count.dtype,np.integer) and (count>=0).all()
 diag={'mean_probability_L1':float(abs((count/libs[:,None]).mean(0)-wp).sum()),'bulk_probability_L1':float(abs(count.sum(0)/libs.sum()-wb).sum()),'conditional_delta_rms':float(np.sqrt(np.mean(delta**2))),'learned_support':int(model['support'][ti])}
 if diag['mean_probability_L1']>.03 or diag['bulk_probability_L1']>.031:raise ValueError('V22 count moments outside tolerance')
 return sparse.csr_matrix(count),diag
