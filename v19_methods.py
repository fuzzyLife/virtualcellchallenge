"""V19 bounded covariance-conditioned response transfer in a training-only program basis."""
from pathlib import Path
import gc,json
import numpy as np
from sklearn.utils.extmath import randomized_svd
import build_submission as b
from fit_v7_response import raw_inputs
from v8_methods import gated
ROOT=b.ROOT;OUT=ROOT/'experiments/v19_covariance_forces';A=OUT/'analysis'
FAMILIES=['H1','HCT116','HEK293T'];NAMES=['v16','covariance','diagonal','shuffled','amplitude']

def lognorm(x,depth):
 y=x.astype(np.float64).multiply((1e4/np.asarray(depth))[:,None]).tocsr();y.data=np.log1p(y.data);return y

def control_stats(x,meta,ids):
 ids=np.asarray(ids);labels=meta.gene_target.to_numpy(str);assert np.all(labels[ids]=='Non-Targeting');batches=meta['sample'].to_numpy(str);names=sorted(set(batches[ids]));bid={v:i for i,v in enumerate(names)};index=np.array([bid[v] for v in batches[ids]])
 depth=meta.total_counts.to_numpy(float);z=np.log1p(depth[ids]);y=lognorm(x[ids],depth[ids]);means=np.array([np.asarray(y[index==i].mean(0)).ravel() for i in range(len(names))]);centers=np.array([z[index==i].mean() for i in range(len(names))]);zd=z-centers[index];counts=np.bincount(index,minlength=len(names));cov=np.asarray(y.T@zd).ravel()/len(ids);vd=max(float(np.mean(zd**2)),1e-12);beta=cov/vd
 variance=np.asarray(y.power(2).mean(0)).ravel()-(means**2*counts[:,None]).sum(0)/len(ids)-cov**2/vd;variance=np.maximum(variance,0)
 return {'y':y,'means':means,'depth_centers':centers,'beta':beta,'variance':variance,'names':names,'index':index,'depth_residual':zd,'control_rows':ids}

def residual_matrix(stat,positions):return stat['y'][:,positions].toarray()-stat['means'][:,positions][stat['index']]-stat['depth_residual'][:,None]*stat['beta'][positions]

def regularized_covariance(scores):
 cov=np.atleast_2d(np.cov(scores,rowvar=False,bias=True));d=np.diag(cov);out=.8*cov+.2*np.diag(d)+.05*max(float(d.mean()),1e-8)*np.eye(len(d));assert np.linalg.eigvalsh(out).min()>0;return out

def bounded_operator(query,source,kind='covariance'):
 q=np.asarray(query);s=np.asarray(source)
 if kind=='diagonal':q=np.diag(np.diag(q));s=np.diag(np.diag(s))
 elif kind=='shuffled':
  perm=np.random.default_rng(19019).permutation(len(q));q=q[np.ix_(perm,perm)]
 elif kind!='covariance':raise ValueError(kind)
 operator=np.linalg.solve(s.T,q.T).T;u,values,vh=np.linalg.svd(operator,full_matrices=False);return (u*np.clip(values,.5,2))@vh

def fit_model(excluded,data,plan):
 training=[f for f in FAMILIES if f not in excluded];assert len(training)>=2;stats={};measures=[];source_records=[]
 contexts={c['family']:c for c in plan['contexts']}
 for family in training:
  x,meta,measured,source=raw_inputs(family,data['genes']);stats[family]=control_stats(x,meta,contexts[family]['left']);measures.append(measured);source_records.append({'family':family,'left_NTC_rows':len(contexts[family]['left']),'source':source});del x;gc.collect()
 variance=np.mean([stats[f]['variance'] for f in training],axis=0);common=np.logical_and.reduce(measures)&~np.isin(data['genes'],data['targets'])&(variance>1e-4);possible=np.flatnonzero(common);positions=possible[np.argsort(variance[possible],kind='stable')[-min(2048,len(possible)):]];assert len(positions)>=256
 scale=np.sqrt(variance[positions]);parts=[residual_matrix(stats[f],positions)/scale/np.sqrt(len(stats[f]['control_rows'])) for f in training];matrix=np.concatenate(parts);_,singular,basis=randomized_svd(matrix,n_components=16,n_iter=3,random_state=19019);basis=basis.T;del matrix,parts
 covariances=[];responses=[];sizes=[]
 for family in training:
  stat=stats[family];control=residual_matrix(stat,positions)/scale@basis;covariances.append(regularized_covariance(control));x,meta,_,_=raw_inputs(family,data['genes']);names={v:i for i,v in enumerate(stat['names'])};batches=meta['sample'].to_numpy(str);depth=meta.total_counts.to_numpy(float);index=np.array([names.get(v,-1) for v in batches]);scores=np.zeros((x.shape[0],16));projection=basis/scale[:,None];mean_projection=stat['means'][:,positions]@projection;beta_projection=stat['beta'][positions]@projection
  for lo in range(0,x.shape[0],2000):
   hi=min(lo+2000,x.shape[0]);block=lognorm(x[lo:hi,positions],depth[lo:hi]);valid=index[lo:hi]>=0;bi=np.maximum(index[lo:hi],0);scores[lo:hi]=block@projection-mean_projection[bi]-(np.log1p(depth[lo:hi])-stat['depth_centers'][bi])[:,None]*beta_projection;scores[lo:hi][~valid]=0
  labels=meta.gene_target.to_numpy(str);mean=np.zeros((300,16));count=np.zeros(300,int)
  for ti,target in enumerate(data['targets']):
   rows=np.flatnonzero((labels==target)&(index>=0));count[ti]=len(rows)
   if len(rows)>=40:mean[ti]=scores[rows].mean(0)
  responses.append(mean);sizes.append(count);del x,scores;gc.collect();print('V19 source programs',excluded,family,'supported',int((count>=40).sum()),flush=True)
 return {'positions':positions,'genes':data['genes'][positions],'scale':scale,'basis':basis,'singular_values':singular,'source_covariances':np.asarray(covariances),'source_shifts':np.asarray(responses),'source_sizes':np.asarray(sizes),'training_families':np.asarray(training,str),'targets':data['targets'],'source_records':source_records}

def query_corrections(model,x,meta,ids,measured):
 positions=model['positions'];coverage=float(measured[positions].mean());result={name:np.zeros((300,x.shape[1]),np.float32) for name in ['covariance','diagonal','shuffled']}
 if coverage<.8:return result,{'coverage':coverage,'abstained':True}
 stat=control_stats(x,meta,ids);scores=residual_matrix(stat,positions)/model['scale']@model['basis'];cov=regularized_covariance(scores)
 y=stat['y'][:,positions].copy();y.data=np.expm1(y.data);y.data=y.data/(1+y.data);derivative=np.asarray(y.mean(0)).ravel()*np.log(2);conversion=model['scale']/np.maximum(derivative,.1);sizes=model['source_sizes'];weights=(sizes>=40)*sizes/(sizes+100);denom=weights.sum(0)
 for kind in result:
  delta=np.zeros((300,16))
  for si in range(len(sizes)):
   op=bounded_operator(cov,model['source_covariances'][si],kind);change=model['source_shifts'][si]@(op-np.eye(16)).T;delta+=weights[si,:,None]*change
  delta=np.divide(delta,denom[:,None],out=np.zeros_like(delta),where=denom[:,None]>0);raw=(delta@model['basis'].T)*conversion;raw[:,~measured[positions]]=0;raw[abs(raw)<1e-8]=0;result[kind][:,positions]=(.25*np.clip(raw,-1,1)).astype(np.float32)
 return result,{'coverage':coverage,'abstained':False,'query_covariance':cov.tolist(),'source_supported_targets':int((denom>0).sum()),'correction_nonzero_targets':int(np.count_nonzero(abs(result['covariance']).sum(1))),'maximum_absolute_correction':float(abs(result['covariance']).max())}

def effects_for_target(base,corrections,ti,query,own):
 effects={'v16':np.asarray(base,np.float32).copy()}
 for kind in ['covariance','diagonal','shuffled']:
  correction=gated(corrections[kind][ti].copy(),np.nan_to_num(query),own);effects[kind]=np.clip(base+correction,-3,3).astype(np.float32);effects[kind][own]=0
 mask=np.arange(len(base))!=own;ratio=float(abs(effects['covariance'][mask]).sum()/max(abs(base[mask]).sum(),1e-12));effects['amplitude']=(base*ratio).astype(np.float32);effects['amplitude'][own]=base[own];return effects

def load_plan():
 plan=json.loads((A/'plan.json').read_text())
 for name,h in plan['dependencies'].items():assert b.sha256(ROOT/name)==h,'Frozen V19 dependency changed: '+name
 return plan

def load_model(family):
 f=A/'models'/f'{family}.npz';rec=json.loads(f.with_suffix('.json').read_text());assert rec['sha256']==b.sha256(f) and rec['plan_sha256']==b.sha256(A/'plan.json')
 with np.load(f,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
