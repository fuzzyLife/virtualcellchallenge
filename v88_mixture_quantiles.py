"""V88: converged covariance-mixture priors and mixture-quantile response projection."""
from pathlib import Path
from functools import lru_cache
import gc
import numpy as np
from scipy.special import logsumexp,ndtr
import v86_joint_posterior as old
from v57_pie import read,write,sha
from v63_identity_guard import guard
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v88_mixture_quantiles/analysis';ARMS=('exact125','shared125','exact125_shuffled')
COMPONENTS=old.COMPONENTS
data=old.data

def check():
 p=read(A/'plan.json')
 for n,h in p['code_sha256'].items():assert sha(ROOT/n)==h,n
 assert sha(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')==p['source_cache_sha256']
 return p

def gaussian(y,var,cov):
 """Unrolled small Cholesky solves; avoid millions of tiny generic LAPACK calls."""
 y=np.asarray(y,float);var=np.asarray(var,float);n,k=y.shape;c=np.broadcast_to(cov,(n,k,k)).copy();ix=np.arange(k);c[:,ix,ix]+=var;l=np.zeros_like(c)
 for i in range(k):
  for j in range(i+1):
   z=c[:,i,j]-np.sum(l[:,i,:j]*l[:,j,:j],axis=1)
   if i==j:assert (z>0).all();l[:,i,j]=np.sqrt(z)
   else:l[:,i,j]=z/l[:,j,j]
 z=np.empty_like(y)
 for i in range(k):z[:,i]=(y[:,i]-np.sum(l[:,i,:i]*z[:,:i],axis=1))/l[:,i,i]
 alpha=np.empty_like(y)
 for i in range(k-1,-1,-1):alpha[:,i]=(z[:,i]-np.sum(l[:,i+1:,i]*alpha[:,i+1:],axis=1))/l[:,i,i]
 inv_l=np.zeros_like(c)
 for i in range(k):
  for j in range(i+1):inv_l[:,i,j]=((1. if i==j else 0.)-np.sum(l[:,i,:i]*inv_l[:,:i,j],axis=1))/l[:,i,i]
 inv=np.einsum('nji,njk->nik',inv_l,inv_l);ld=2*np.log(l[:,ix,ix]).sum(1);ll=-.5*(ld+np.sum(z*z,axis=1)+k*np.log(2*np.pi));return ll,inv,alpha

def likelihood(y,var,valid,sim,query_similarity=None):
 n=len(y);k=len(COMPONENTS);ll=np.zeros((n,k));means=np.zeros_like(ll);variances=np.zeros_like(ll);bits=(valid.astype(np.uint16)*(1<<np.arange(y.shape[1]))).sum(1)
 for bit in np.unique(bits):
  rows=np.flatnonzero(bits==bit);indices=np.flatnonzero((int(bit)&(1<<np.arange(y.shape[1])))>0)
  if not len(indices):continue
  yy=y[rows][:,indices];vv=var[rows][:,indices]
  for ci,(scale,length) in enumerate(COMPONENTS):
   cov=np.zeros((len(indices),len(indices))) if not scale else scale**2*old.kernel(sim,length)[np.ix_(indices,indices)]
   values,inverse,alpha=gaussian(yy,vv,cov);ll[rows,ci]=values
   if query_similarity is not None and scale:
    cross=scale**2*.8*(np.ones(len(indices)) if length==999 else np.exp(-np.maximum(1-query_similarity[indices],0)/length));means[rows,ci]=alpha@cross;variances[rows,ci]=np.maximum(scale**2-np.einsum('i,nij,j->n',cross,inverse,cross),0)
 return ll,means,variances

def objective(a,w):return float(np.log(a@w).sum()+.5*np.log(w).sum())
def solve_weights(ll):
 """Strictly concave positive-simplex MAP mixture fit with a KKT stopping rule."""
 a=np.exp(ll-np.max(ll,axis=1,keepdims=True));n,k=a.shape;w=np.full(k,1/k);history=[];lam=n+.5*k
 for it in range(150):
  t=a@w;g=a.T@(1/t)+.5/w;res=float(np.max(abs(g-lam))/lam);obj=objective(a,w);history.append(dict(iteration=it,objective_shifted=obj,kkt_relative_residual=res))
  if res<1e-8:break
  h=(a/t[:,None]).T@(a/t[:,None])+np.diag(.5/w**2);system=np.zeros((k+1,k+1));system[:k,:k]=h;system[:k,k]=1;system[k,:k]=1;step=np.linalg.solve(system,np.r_[g,0])[:k];rate=min(1.,.99*np.min(-w[step<0]/step[step<0]) if np.any(step<0) else 1.);slope=float(g@step)
  for backtrack in range(70):
   proposal=w+rate*step
   if (proposal>0).all() and objective(a,proposal)>=obj+1e-4*rate*slope-1e-9:break
   rate*=.5
  else:raise RuntimeError('Newton line search failed')
  w=proposal;w/=w.sum()
 else:raise RuntimeError('Mixture weight optimization did not converge')
 assert (w>0).all() and abs(w.sum()-1)<1e-12 and history[-1]['kkt_relative_residual']<1e-8
 return w,history

def probabilities(ll,w):
 z=ll+np.log(np.maximum(w,1e-300))[None,:];return np.exp(z-logsumexp(z,axis=1)[:,None])
def cdf(x,mu,var,prob,left=False):
 sd=np.sqrt(var);z=np.divide(x[:,None]-mu,sd,out=np.zeros_like(mu),where=sd>0);value=ndtr(z);atom=sd==0;value[atom]=np.broadcast_to((x[:,None]>mu) if left else (x[:,None]>=mu),mu.shape)[atom];return np.sum(prob*value,axis=1)

def quantile(mu,var,prob,q):
 """Numerical generalized-inverse CDF, explicitly retaining the null point mass."""
 zero=np.zeros(len(mu));below=cdf(zero,mu,var,prob,True);at=cdf(zero,mu,var,prob);answer=np.zeros(len(mu));active=~((below<q)&(at>=q));ids=np.flatnonzero(active)
 if len(ids):
  mm=mu[ids];vv=var[ids];pp=prob[ids];sd=np.sqrt(vv);lo=np.min(mm-12*sd,axis=1);hi=np.max(mm+12*sd,axis=1)
  for _ in range(32):
   mid=(lo+hi)*.5;less=cdf(mid,mm,vv,pp)<q;lo=np.where(less,mid,lo);hi=np.where(less,hi,mid)
  answer[ids]=hi
 return answer

def project(base,mu,var,ll,w):
 prob=probabilities(ll,w);below=cdf(base,mu,var,prob,True);at=cdf(base,mu,var,prob);answer=base.copy();low=at<.1;high=below>.9
 if low.any():answer[low]=quantile(mu[low],var[low],prob[low],.1)
 if high.any():answer[high]=quantile(mu[high],var[high],prob[high],.9)
 return answer

def fit(outer,shuffle=False):
 check();b=old.bank(outer);y,v,ok=old.observations(b,shuffle);eligible=np.flatnonzero(ok.sum(1)>=2);ids=np.random.default_rng(86086).choice(eligible,min(20000,len(eligible)),replace=False);ll,_,_=likelihood(y[ids],v[ids],ok[ids],b['similarity']);w,h=solve_weights(ll);shared_ix=np.array([i for i,(s,l) in enumerate(COMPONENTS) if s==0 or l==999]);sw,sh=solve_weights(ll[:,shared_ix]);shared=np.zeros(len(w));shared[shared_ix]=sw
 folder=A/'models';folder.mkdir(exist_ok=True);p=folder/(outer+('_shuffled' if shuffle else '')+'.npz');assert not p.exists();np.savez_compressed(p,weights=w,shared_weights=shared,sample_indices=ids,axis=b['axis'],rank_features=b['rank_features'],similarity=b['similarity'],source_names=b['source_names']);comparison=None
 if outer!='production':
  op=old.A/'models'/p.name;om=read(op.with_suffix('.json'));assert sha(op)==om['sha256']
  with np.load(op,allow_pickle=False) as z:np.testing.assert_array_equal(z['sample_indices'],ids);old_w=z['weights'].copy()
  a=np.exp(ll-ll.max(1,keepdims=True));comparison=dict(V86_shifted_objective=objective(a,old_w),converged_shifted_objective=objective(a,w),maximum_weight_difference=float(abs(w-old_w).max()))
 write(p.with_suffix('.json'),dict(sha256=sha(p),outer=outer,shuffled=shuffle,source_names=b['source_names'].tolist(),weights=w.tolist(),shared_weights=shared.tolist(),optimization=h,shared_optimization=sh,comparison_to_V86=comparison,plan_sha256=sha(A/'plan.json')));print('V88 fit',outer,shuffle,'steps',len(h),'KKT',h[-1]['kkt_relative_residual'],flush=True);del y,v,ok,ll;gc.collect()

@lru_cache(7)
def model(outer,shuffle):
 p=A/'models'/(outer+('_shuffled' if shuffle else '')+'.npz');assert sha(p)==read(p.with_suffix('.json'))['sha256']
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def response(query,base,outer,shuffle=False,controls=True):
 b=old.bank(outer);net=model(outer,shuffle);axis=net['axis'];observed=np.isfinite(query[axis]);assert observed.mean()>=.9;feature=np.zeros(len(axis));r=__import__('scipy').stats.rankdata(query[axis][observed]);r-=r.mean();feature[observed]=r/np.linalg.norm(r);sim=net['rank_features']@feature;y,v,ok=old.observations(b,shuffle);rows=np.flatnonzero(ok.any(1));flat=base.ravel();dest=flat.copy();shared=flat.copy()
 for start in range(0,len(rows),30000):
  ii=rows[start:start+30000];ll,mu,vv=likelihood(y[ii],v[ii],ok[ii],net['similarity'],sim);dest[ii]=project(flat[ii],mu,vv,ll,net['weights'])
  if controls and not shuffle:shared[ii]=project(flat[ii],mu,vv,ll,net['shared_weights'])
  if start%990000==0:print('V88 posterior progress',outer,'shuffled',shuffle,start,'/',len(rows),flush=True)
 return dict(exact=dest.reshape(base.shape),shared=shared.reshape(base.shape),supported=ok.any(1).reshape(base.shape)),dict(source_names=net['source_names'].tolist(),rank_inner_products=sim.tolist(),measured_rank_fraction=float(observed.mean()))

def predict(query,base,bulk,outer,controls=True):
 check();d=old.data();own=np.array([list(d['genes']).index(t) for t in d['targets']]);panel=np.isin(d['genes'],d['targets']);assay=np.log1p(.05*np.nan_to_num(query)[None,:]*np.exp2(base))-np.log1p(.05*np.nan_to_num(query)[None,:]);out={};trace={}
 for fake in ([False,True] if controls else [False]):
  post,tr=response(query,assay,outer,fake,controls)
  for typ in (['exact'] if fake or not controls else ['exact','shared']):
   name=typ+'125'+('_shuffled' if fake else '');desired=old.invert_assay(post[typ],query);mask=post['supported']&~panel[None,:]&np.isfinite(query)[None,:]&(query[None,:]>0);correction=np.where(mask,.125*np.clip(desired-base,-.5,.5),0);proposal=np.clip(base+correction,-3,3).astype(np.float32);proposal[np.arange(300),own]=0;value,diag=guard(base,proposal,bulk,own,panel);inactive=~np.any(correction!=0,axis=1);value[inactive]=base[inactive];out[name]=value;trace[name]=dict(diag,routing=tr,changed_targets=int(np.any(value!=base,axis=1).sum()))
 return out,trace
