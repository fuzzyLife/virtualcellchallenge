"""V25 intervention-aligned response geometry; experimental hypothesis, not a physical law."""
from pathlib import Path
import json,gc
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy import stats
from scipy import sparse
from sklearn.utils.extmath import randomized_svd
import build_submission as b
from fit_v7_response import read_cache,raw_inputs
from v6_methods import moments
from v16_methods import source_arrays,emit as old_emit
from v24_methods import make_teacher,predict as mean_predict
from v19_methods import control_stats,residual_matrix
from v8_methods import gated
from v22_methods import write,digest,logcpm
from analyze_target_rank_correlations import pairwise_spearman
ROOT=b.ROOT;OUT=ROOT/'experiments/v25_cross_domain';A=OUT/'analysis'
FAMILIES=['H1','HCT116','HEK293T','CD4'];NAMES=['v24','force25','force50','bures25','diagonal25','shuffled25','shuffled50','amplitude25','source_program25'];ELIGIBLE=['force25','force50','bures25']
CONTEXTS=[('H1',4),('HCT116',5),('HEK293T',6),('CD4_Rest_D2',1),('CD4_Stim48_D2',3)]
def load_plan():
 p=json.loads((A/'plan.json').read_text())
 for n,h in p['dependencies'].items():assert b.sha256(ROOT/n)==h,'Frozen V25 code changed: '+n
 for n,h in p['source_inputs'].items():assert b.sha256(ROOT/n)==h
 return p

def raw(family,genes):
 if family in FAMILIES[:3]:return raw_inputs(family,genes)
 folders={'CD4_Rest_D1':'experiments/v22_state_response/analysis/cd4_D1_selected','CD4_Rest_D2':'experiments/v22_state_response/analysis/cd4_selected','CD4_Stim48_D2':'experiments/v24_cd4_response/analysis/holdout','CD4_D3_Rest':'experiments/v25_cross_domain/analysis/holdout'};folder=ROOT/folders[family];rec=json.loads((folder/'extraction.json').read_text());x=sparse.load_npz(folder/'counts.npz');meta=pd.read_csv(folder/'cells.csv');axis=pd.read_csv(folder/'genes.csv');assert b.sha256(folder/'counts.npz')==rec['counts_sha256'] and b.sha256(folder/'cells.csv')==rec['cells_sha256'];np.testing.assert_array_equal(axis.gene_name.to_numpy(str),genes);meta['sample']=meta.lane_id.astype(str);return x,meta,axis.measured.to_numpy(bool),rec

def covariance(z):
 c=np.atleast_2d(np.cov(z,rowvar=False,bias=True));d=np.maximum(np.diag(c),1e-8);return .75*c+.25*np.diag(d)+.1*d.mean()*np.eye(len(d))
def power(c,a):
 e,u=np.linalg.eigh((c+c.T)/2);assert e.min()>0;return (u*(e**a))@u.T

def operator(q,s,kind):
 if kind=='diagonal':q=np.diag(np.diag(q));s=np.diag(np.diag(s))
 if kind=='shuffled':
  perm=np.random.default_rng(25025).permutation(len(q));q=q[np.ix_(perm,perm)]
 if kind=='bures':
  half=power(s,.5);inverse=power(s,-.5);op=inverse@power(half@q@half,.5)@inverse
 else:op=np.linalg.solve(s.T,q.T).T
 u,e,v=np.linalg.svd(op,full_matrices=False);return (u*np.clip(e,2/3,1.5))@v

def force_confidence(a,c):
 return (np.sign(a)==np.sign(c))*np.divide(2*np.minimum(abs(a),abs(c)),abs(a)+abs(c),out=np.zeros_like(a),where=abs(a)+abs(c)>1e-10)

def fit_model(excluded):
 plan=load_plan();data=read_cache();genes=data['genes'];targets=data['targets'];sources=[];stat=[];records=[];means=[];valids=[];reliability=[];profiles=[];measured=[]
 for family,ci in CONTEXTS:
  f='CD4' if family.startswith('CD4') else family
  if f in excluded:continue
  x,m,meas,rec=raw(family,genes);ids=np.flatnonzero(m.gene_target=='Non-Targeting');st=control_stats(x,m,ids);stat.append(st);measured.append(meas);records.append({'context':family,'family':f,'control_rows':ids.tolist(),'source':rec});mu=np.zeros((300,len(genes)));den=np.zeros_like(mu);semean=np.zeros_like(mu)
  for ri in np.flatnonzero(data['row_contexts']==ci):
   ti=int(data['row_targets'][ri]);mask=data['valid'][ri];own=int(np.flatnonzero(genes==targets[ti])[0]);kd=float(data['knockdown'][ri]);kd=.9 if not np.isfinite(data['profiles'][ci,own]) or data['profiles'][ci,own]<=5 else kd;gain=float(np.clip(.9/max(kd,.35),.8,1.8));var=np.where(mask,np.maximum((data['errors'][ri]*gain)**2,.0025),1e6);w=np.where(mask,np.clip(.04/var,.1,4),0);mu[ti]+=w*np.where(mask,np.clip(data['truth'][ri]*gain,-3,3),0);den[ti]+=w;semean[ti]+=w*var
  mu=np.divide(mu,den,out=np.zeros_like(mu),where=den>0);var=np.divide(semean,den,out=np.ones_like(mu),where=den>0);valid=den>0;profile=data['profiles'][ci];derivative=np.log(2)*np.nan_to_num(profile)*.01/(1+np.nan_to_num(profile)*.01);means.append(mu*derivative);valids.append(valid);reliability.append(np.clip(.04/np.maximum(var,.0025),.1,4));profiles.append(profile);sources.append(family);del x;gc.collect();print('V25 training controls read',excluded,family,flush=True)
 variance=np.mean([s['variance'] for s in stat],axis=0);common=np.logical_and.reduce(measured)&~np.isin(genes,targets)&(variance>1e-4)&np.isfinite(profiles).all(0);common&=np.logical_and.reduce([v.any(0) for v in valids]);scale_all=np.sqrt(np.maximum(variance,1e-4));energy=np.mean([np.mean(y*y*v,axis=0) for y,v in zip(means,valids)],axis=0)/(variance+.01);possible=np.flatnonzero(common);positions=possible[np.argsort(energy[possible],kind='stable')[-min(4096,len(possible)):]];assert len(positions)>=512;scale=scale_all[positions];parts=[]
 for y,v in zip(means,valids):
  coverage=v[:,positions].mean(1);keep=coverage>=.9;parts.append(y[keep][:,positions]/scale/np.sqrt(max(keep.sum(),1)))
 matrix=np.concatenate(parts);_,singular,basis=randomized_svd(matrix,n_components=64,n_iter=4,random_state=25025);basis=basis.T;covs=[];shifts=[];weights=[]
 for i,s in enumerate(stat):
  score=residual_matrix(s,positions)/scale@basis;covs.append(covariance(score));coverage=valids[i][:,positions].mean(1);delta=means[i][:,positions]/scale@basis;delta[coverage<.9]=0;shifts.append(delta);precision=np.mean(reliability[i][:,positions],axis=1);family_weight={'H1':2.,'HCT116':1.,'HEK293T':1.,'CD4_Rest_D2':.125,'CD4_Stim48_D2':.125}[sources[i]];weights.append(family_weight*precision*coverage*(coverage>=.9))
 result={'genes':genes.astype(str),'targets':targets.astype(str),'positions':positions,'scale':scale,'basis':basis,'singular_values':singular,'source_covariances':np.array(covs),'source_shifts':np.array(shifts),'source_weights':np.array(weights),'source_profiles':np.array(profiles),'source_contexts':np.array(sources,str)};folder=A/'models';folder.mkdir(exist_ok=True);file=folder/(('_'.join(excluded) or 'all')+'.npz');assert not file.exists();np.savez_compressed(file,**result);write(file.with_suffix('.json'),{'sha256':b.sha256(file),'plan_sha256':b.sha256(A/'plan.json'),'excluded_families':excluded,'training_contexts':sources,'source_records':records,'basis':'64 training intervention-aligned directions on up to4096 shared nonpanel output genes. Missing source-response coordinates zero in projection only when>=90% measured, otherwise abstain; not biological zero votes. LogCP10k derivative linearization; empirical Gaussian covariance proxy, no asserted universal law.'});print('V25 response geometry fitted',excluded,len(positions),flush=True)

def load_model(family):
 key='CD4' if family.startswith('CD4') else family;key='all' if family=='all' else key;file=A/'models'/(key+'.npz');rec=json.loads(file.with_suffix('.json').read_text());assert b.sha256(file)==rec['sha256'] and b.sha256(A/'plan.json')==rec['plan_sha256']
 with np.load(file,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def baseline(data,query,excluded=(),production=False,parent=None):
 own=np.array([np.flatnonzero(data['genes']==t)[0] for t in data['targets']])
 if parent is not None:return parent.copy()
 if 'CD4' in excluded:
  bank=source_arrays(data,list(excluded));return np.asarray([gated(v.copy(),np.nan_to_num(query),int(o)) for v,o in zip(bank['base'],own)],np.float32)
 teacher=make_teacher(True);result,_=mean_predict(data,teacher,query,excluded,production=True);return result['cdf50_de']

def predict(data,model,x,meta,ids,measured,excluded=(),parent=None):
 query=moments(x[ids])[1]*1e6;query[~measured]=np.nan;base=baseline(data,query,excluded,parent=parent);pos=model['positions'];result={n:base.copy() for n in NAMES};coverage=float(measured[pos].mean())
 if coverage<.9:return result,{'abstained':True,'coverage':coverage}
 st=control_stats(x,meta,ids);z=residual_matrix(st,pos)/model['scale']@model['basis'];half0=[];half1=[]
 for batch in range(len(st['names'])):
  rows=np.flatnonzero(st['index']==batch);half0.extend(rows[::2]);half1.extend(rows[1::2])
 cq=covariance(z);halves=[covariance(z[half0]),covariance(z[half1])];axis=~np.isin(data['genes'],data['targets']);rho=np.array([pairwise_spearman(query[axis],p[axis],min_pairs=100)[0] for p in model['source_profiles']]);rho=np.nan_to_num(rho,nan=-1);weights=model['source_weights']*np.exp(2*(rho-rho.max()))[:,None];den=weights.sum(0);own=np.array([np.flatnonzero(data['genes']==t)[0] for t in data['targets']]);der=np.maximum(np.log(2)*np.nan_to_num(query[pos])*.01/(1+np.nan_to_num(query[pos])*.01),.1);base_score=(base[:,pos]*der/model['scale'])@model['basis'];inverse=np.linalg.inv(cq);base_norm=np.sqrt(np.maximum(np.einsum('ti,ij,tj->t',base_score,inverse,base_score),0));corrections={};diagnostic={'abstained':False,'coverage':coverage,'source_rho':rho.tolist(),'query_covariance_eigenvalues':np.linalg.eigvalsh(cq).tolist()}
 for kind in ['force','bures','diagonal','shuffled']:
  delta=[]
  for q in halves:
   total=np.zeros((300,64))
   for si,s in enumerate(model['source_covariances']):total+=weights[si,:,None]*(model['source_shifts'][si]@(operator(q,s,kind)-np.eye(64)).T)
   delta.append(np.divide(total,den[:,None],out=np.zeros_like(total),where=den[:,None]>0))
  change=.5*(delta[0]+delta[1])*force_confidence(delta[0],delta[1]);norm=np.sqrt(np.maximum(np.einsum('ti,ij,tj->t',change,inverse,change),0));bound=np.minimum(1,np.divide(.5*base_norm,np.maximum(norm,1e-12)));change*=bound[:,None];raw_delta=(change@model['basis'].T)*model['scale']/der;raw_delta[abs(raw_delta)<1e-8]=0;corrections[kind]=raw_delta;diagnostic[kind+'_mean_split_agreement']=float(force_confidence(delta[0],delta[1]).mean())
 for name,kind,alpha in [('force25','force',.25),('force50','force',.5),('bures25','bures',.25),('diagonal25','diagonal',.25),('shuffled25','shuffled',.25),('shuffled50','shuffled',.5)]:
  extra=np.zeros_like(base);extra[:,pos]=alpha*np.clip(corrections[kind],-1,1);extra=np.asarray([gated(v,np.nan_to_num(query),int(o)) for v,o in zip(extra,own)],np.float32);result[name]=np.clip(base+extra,-3,3).astype(np.float32)
 average=np.divide((weights[:,:,None]*model['source_shifts']).sum(0),den[:,None],out=np.zeros((300,64)),where=den[:,None]>0);projected=(average@model['basis'].T)*model['scale']/der;extra=np.zeros_like(base);extra[:,pos]=.25*np.clip(projected-base[:,pos],-1,1);extra=np.asarray([gated(v,np.nan_to_num(query),int(o)) for v,o in zip(extra,own)],np.float32);result['source_program25']=np.clip(base+extra,-3,3).astype(np.float32);ratio=np.divide(abs(result['force25']).sum(1),np.maximum(abs(base).sum(1),1e-12));result['amplitude25']=np.clip(base*ratio[:,None],-3,3).astype(np.float32)
 for v in result.values():v[np.arange(300),own]=0
 return result,diagnostic

def emit(template,effect,own,seed,prepared,atlas):return old_emit(template,np.arange(template.shape[0]),effect,own,seed,'residual10',prepared,atlas)
