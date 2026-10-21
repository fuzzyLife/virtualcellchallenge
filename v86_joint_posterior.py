"""V86 covariance-mixture Gaussian-process response transfer, inspired by multivariate shrinkage."""
from pathlib import Path
from functools import lru_cache
import gc
import numpy as np
from scipy import stats
from scipy.special import logsumexp
from v57_pie import read,write,sha
from v15_ensemble_methods import cube,profiles as source_profiles
from v7_response_methods import FAMILIES
from v62_conserved_response import data,invert_assay
from v63_identity_guard import guard

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'experiments/v86_joint_posterior';A=OUT/'analysis'
ARMS=('interval25','posterior25','shared25','interval25_shuffled')
GROUP={'K562':'Replogle','CD4':'CD4','H1':'H1','HCT116':'Orion','HEK293T':'Orion'}
COMPONENTS=[(0.,-1.)]+[(s,l) for s in [.025,.1,.4,1.] for l in [.05,.2,1.,4.,999.]]

def check():
    p=read(A/'plan.json')
    for n,h in p['code_sha256'].items():assert sha(ROOT/n)==h,n
    return p

@lru_cache(3)
def bank(outer):
    d=data();excluded=[str(f) for f in FAMILIES if GROUP[str(f)]==outer];y,v,ok,own=cube(d,excluded);pp=source_profiles(d);keep=np.array([i for i,f in enumerate(FAMILIES) if str(f) not in excluded]);y=y[keep];v=v[keep];ok=ok[keep];pp=pp[keep]
    amount=.05*np.nan_to_num(pp)[:,None,:];derivative=np.log(2)*amount*np.exp2(y)/(1+amount*np.exp2(y));yy=np.log1p(amount*np.exp2(y))-np.log1p(amount);vv=np.maximum(v*derivative**2,.0025)
    nonpanel=~np.isin(d['genes'],d['targets']);ok&=nonpanel[None,None,:];axis=np.flatnonzero(nonpanel&np.isfinite(pp).all(0)&(pp>5).all(0));assert len(axis)>=100
    rank=stats.rankdata(pp[:,axis],axis=1,method='average');rank-=rank.mean(1,keepdims=True);rank/=np.linalg.norm(rank,axis=1,keepdims=True);similarity=rank@rank.T
    return dict(y=yy.astype(np.float32),variance=vv.astype(np.float32),valid=ok,source_profiles=pp,rank_features=rank,axis=axis,similarity=similarity,source_indices=keep,source_names=FAMILIES[keep].astype(str))

def observations(b,shuffle=False):
    y=b['y'];v=b['variance'];ok=b['valid']
    if shuffle:
        n=y.shape[1];order=np.array([np.roll(np.arange(n),37+11*i) for i in range(len(y))]);y=np.array([row[ix] for row,ix in zip(y,order)]);v=np.array([row[ix] for row,ix in zip(v,order)]);ok=np.array([row[ix] for row,ix in zip(ok,order)])
    return y.reshape(len(y),-1).T,v.reshape(len(y),-1).T,ok.reshape(len(y),-1).T

def kernel(similarity,length):
    base=np.ones_like(similarity) if length==999 else np.exp(-np.maximum(1-similarity,0)/length)
    return .8*base+.2*np.eye(len(base))

def component_likelihood(y,var,cov):
    c=np.broadcast_to(cov,(len(y),)+cov.shape).copy();ix=np.arange(y.shape[1]);c[:,ix,ix]+=var
    sign,ld=np.linalg.slogdet(c);assert (sign>0).all();inv=np.linalg.inv(c);alpha=np.einsum('nij,nj->ni',inv,y)
    ll=-.5*(ld+np.einsum('ni,ni->n',y,alpha)+y.shape[1]*np.log(2*np.pi))
    return ll,inv,alpha

def likelihood(y,var,valid,sim,query_similarity=None):
    n=len(y);k=len(COMPONENTS);ll=np.zeros((n,k));means=np.zeros_like(ll);variances=np.zeros_like(ll);bits=(valid.astype(np.uint16)*(1<<np.arange(y.shape[1]))).sum(1)
    for bit in np.unique(bits):
        rows=np.flatnonzero(bits==bit);indices=np.flatnonzero((int(bit)&(1<<np.arange(y.shape[1])))>0)
        if not len(indices):continue
        yy=y[rows][:,indices];vv=var[rows][:,indices]
        for ci,(scale,length) in enumerate(COMPONENTS):
            cov=np.zeros((len(indices),len(indices))) if not scale else scale**2*kernel(sim,length)[np.ix_(indices,indices)]
            values,inverse,alpha=component_likelihood(yy,vv,cov);ll[rows,ci]=values
            if query_similarity is not None and scale:
                cross=scale**2*.8*(np.ones(len(indices)) if length==999 else np.exp(-np.maximum(1-query_similarity[indices],0)/length))
                means[rows,ci]=alpha@cross;variances[rows,ci]=np.maximum(scale**2-np.einsum('i,nij,j->n',cross,inverse,cross),0)
    return ll,means,variances

def em(ll):
    w=np.full(ll.shape[1],1/ll.shape[1]);history=[]
    for iteration in range(100):
        total=ll+np.log(w)[None,:];norm=logsumexp(total,axis=1);r=np.exp(total-norm[:,None]);new=(r.sum(0)+.5)/(len(ll)+.5*ll.shape[1]);change=float(abs(new-w).max());w=new;history.append(dict(iteration=iteration,mean_log_likelihood=float(norm.mean()),maximum_weight_change=change))
        if change<1e-7:break
    return w,history

def fit(outer,shuffle=False):
    check();b=bank(outer);y,v,ok=observations(b,shuffle);eligible=np.flatnonzero(ok.sum(1)>=2);assert len(eligible)>=1000
    ids=np.random.default_rng(86086).choice(eligible,min(20000,len(eligible)),replace=False);ll,_,_=likelihood(y[ids],v[ids],ok[ids],b['similarity']);w,h=em(ll);indices=np.array([i for i,(s,l) in enumerate(COMPONENTS) if s==0 or l==999]);sw,sh=em(ll[:,indices]);shared=np.zeros(len(w));shared[indices]=sw
    dest=A/'models';dest.mkdir(exist_ok=True);p=dest/(outer+('_shuffled' if shuffle else '')+'.npz');assert not p.exists();np.savez_compressed(p,weights=w,shared_weights=shared,sample_indices=ids,components=np.asarray(COMPONENTS),similarity=b['similarity'],axis=b['axis'],rank_features=b['rank_features'],source_names=b['source_names'])
    write(p.with_suffix('.json'),dict(sha256=sha(p),outer=outer,shuffled_training_and_conditioning=shuffle,source_names=b['source_names'].tolist(),observed_coordinates=len(eligible),fitting_coordinates=len(ids),weights=w.tolist(),shared_weights=shared.tolist(),em_history=h,shared_em_history=sh,plan_sha256=sha(A/'plan.json'),source_cache_sha256=sha(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),script_sha256=sha(__file__),scope='Correlated gene/target coordinates are used as working-likelihood training units, not independent biological samples. Whole recipient study excluded.'))
    print('V86 model fitted',outer,'shuffled',shuffle,'n',len(ids),flush=True);del y,v,ok,ll;gc.collect()

@lru_cache(6)
def model(outer,shuffle):
    p=A/'models'/(outer+('_shuffled' if shuffle else '')+'.npz');assert sha(p)==read(p.with_suffix('.json'))['sha256']
    with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def mixture(ll,means,var,w):
    total=ll+np.log(np.maximum(w,1e-300))[None,:];prob=np.exp(total-logsumexp(total,axis=1)[:,None]);mu=np.sum(prob*means,axis=1);variance=np.maximum(np.sum(prob*(var+means**2),axis=1)-mu**2,0)
    return mu,variance

def posterior(query,outer,shuffle=False):
    b=bank(outer);net=model(outer,shuffle);axis=net['axis'];observed=np.isfinite(query[axis]);assert observed.mean()>=.9,'Insufficient measured query ranks; explicit abstention required'
    feature=np.zeros(len(axis));rank=stats.rankdata(query[axis][observed],method='average');rank-=rank.mean();den=np.linalg.norm(rank);assert den>0;feature[observed]=rank/den;similarity=net['rank_features']@feature
    y,v,ok=observations(b,shuffle);n=len(y);mu=np.zeros(n);var=np.zeros(n);smu=np.zeros(n);svar=np.zeros(n);active=np.flatnonzero(ok.any(1))
    for start in range(0,len(active),40000):
        rows=active[start:start+40000];ll,means,variance=likelihood(y[rows],v[rows],ok[rows],net['similarity'],similarity);mu[rows],var[rows]=mixture(ll,means,variance,net['weights']);smu[rows],svar[rows]=mixture(ll,means,variance,net['shared_weights'])
    shape=b['y'].shape[1:]
    return dict(mean=mu.reshape(shape),variance=var.reshape(shape),shared_mean=smu.reshape(shape),shared_variance=svar.reshape(shape),supported=ok.any(1).reshape(shape)),dict(query_source_rank_inner_products=similarity.tolist(),query_measured_rank_fraction=float(observed.mean()),rank_axis_genes=len(axis),source_names=b['source_names'].tolist())

def predict(query,base,bulk,outer):
    check();d=data();own=np.array([list(d['genes']).index(t) for t in d['targets']]);panel=np.isin(d['genes'],d['targets']);
    if np.isfinite(query[bank(outer)['axis']]).mean()<.9:return {a:base.copy() for a in ARMS},{a:dict(abstained_insufficient_query_measurement=True) for a in ARMS}
    q=np.nan_to_num(query);base_assay=np.log1p(.05*q[None,:]*np.exp2(base))-np.log1p(.05*q[None,:]);real,tr=posterior(query,outer);fake,ft=posterior(query,outer,True);out={};trace={}
    for name,post,interval,shared in [('interval25',real,True,False),('posterior25',real,False,False),('shared25',real,True,True),('interval25_shuffled',fake,True,False)]:
        mu=post['shared_mean'] if shared else post['mean'];var=post['shared_variance'] if shared else post['variance'];dest=np.clip(base_assay,mu-1.2815515655446004*np.sqrt(var),mu+1.2815515655446004*np.sqrt(var)) if interval else mu
        response=invert_assay(dest,query);mask=post['supported']&~panel[None,:]&np.isfinite(query)[None,:]&(query[None,:]>0);correction=np.where(mask,.25*np.clip(response-base,-.5,.5),0);proposal=np.clip(base+correction,-3,3).astype(np.float32);proposal[np.arange(len(own)),own]=0;value,diag=guard(base,proposal,bulk,own,panel);inactive=~np.any(correction!=0,axis=1);value[inactive]=base[inactive];out[name]=value;trace[name]=dict(diag,changed_targets=int(np.any(value!=base,axis=1).sum()),query_routing=ft if name.endswith('shuffled') else tr)
    return out,trace
