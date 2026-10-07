"""Masked cross-context response pooling with per-target empirical Bayes shrinkage."""
import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm
from analyze_target_rank_correlations import pairwise_spearman
FAMILIES=np.array(['K562','CD4','H1','HCT116','HEK293T'])
BASE_WEIGHTS=np.array([2.,.5,2.,1.,1.])
CONFIGS=[{'name':f'{mode}_a{int(amp*100):03d}','mode':mode,'amplitude':amp}
         for mode in ['uniform','rank','adaptive','adaptive_kd'] for amp in [.35,.65,1.]]
def collapse(data,train):
    nt=len(data['targets']);ng=len(data['genes'])
    y=np.zeros((len(FAMILIES),nt,ng),np.float32)
    var=np.ones_like(y)*1e6;valid=np.zeros_like(y,bool);kd=np.zeros((len(FAMILIES),nt))
    for fi,family in enumerate(FAMILIES):
        for ti in range(nt):
            rows=np.flatnonzero(train&(data['row_families']==family)&(data['row_targets']==ti))
            if not len(rows):continue
            masks=data['valid'][rows];den=masks.sum(axis=0);ok=den>0
            mean=np.divide((data['truth'][rows]*masks).sum(axis=0),den,out=np.zeros(ng),where=ok)
            # CD4 conditions are correlated and heterogeneous: no artificial SE/sqrt(3).
            uncertainty=np.divide(((data['errors'][rows]**2+(data['truth'][rows]-mean)**2)*masks).sum(axis=0),
                den,out=np.ones(ng)*1e6,where=ok)
            y[fi,ti]=mean;var[fi,ti]=np.maximum(uncertainty,.05**2);valid[fi,ti]=ok
            kd[fi,ti]=np.mean(data['knockdown'][rows])
    return y,var,valid,kd
def context_weights(data,query_context,train,rank=True):
    weights=BASE_WEIGHTS.copy();rho=[]
    for family in FAMILIES:
        contexts=np.unique(data['row_contexts'][train&(data['row_families']==family)])
        values=[]
        for ci in contexts:
            r,n,status=pairwise_spearman(data['profiles'][query_context],data['profiles'][ci],min_pairs=100)
            if status=='ok':values.append(r)
        rho.append(np.mean(values) if values else np.nan)
    rho=np.asarray(rho)
    if rank:
        available=np.isfinite(rho)
        weights[available]*=np.exp(4*(rho[available]-np.max(rho[available])))
    weights[~np.isfinite(rho)]=0
    return weights,rho
def combine(cube,weights,normalize_kd=False):
    y,var,valid,kd=cube
    y=y.copy()
    if normalize_kd:
        correction=np.clip(.9/np.maximum(kd,.25),.75,1.5)
        y*=correction[:,:,None];var=var*correction[:,:,None]**2
    # Measurement confidence bounded so a single nearly noiseless source cannot dominate.
    confidence=np.clip(.04/np.maximum(var,.0025),.1,4.)
    w=weights[:,None,None]*confidence*valid
    den=w.sum(axis=0);ok=den>0
    wn=np.divide(w,den[None,:],out=np.zeros_like(w),where=den[None,:]>0)
    mean=(wn*y).sum(axis=0)
    measurement=(wn*wn*var).sum(axis=0)
    heterogeneity=(wn*(y-mean[None,:])**2).sum(axis=0)
    # Retain across-source disagreement as transfer uncertainty; never treat masks as zero votes.
    variance=np.maximum(measurement+.25*heterogeneity,.0025)
    mean[~ok]=0;variance[~ok]=1e6
    return mean.astype(np.float32),variance.astype(np.float32),ok
def adaptive_shrink(mean,variance,valid):
    """Target-specific mixture-of-normal prior estimated using training-source responses."""
    out=np.zeros_like(mean);priors=[]
    scales=np.array([0.,.05,.15,.4,1.,2.])
    for ti in range(len(mean)):
        columns=np.flatnonzero(valid[ti])
        if not len(columns):priors.append(np.zeros(len(scales)));continue
        # Deterministic thinning bounds EM cost, independent of held-out responses.
        fit=columns[::max(1,len(columns)//2500)]
        value=mean[ti,fit].astype(float);noise=variance[ti,fit].astype(float)
        total=noise[:,None]+scales[None,:]**2
        likelihood=-.5*(np.log(total)+value[:,None]**2/total)
        mix=np.ones(len(scales))/len(scales)
        for _ in range(30):
            logp=likelihood+np.log(np.maximum(mix,1e-12))
            responsibilities=np.exp(logp-logsumexp(logp,axis=1)[:,None])
            new=(responsibilities.sum(axis=0)+.1)/(len(fit)+.1*len(scales))
            if np.max(np.abs(new-mix))<1e-5:mix=new;break
            mix=new
        value=mean[ti,columns].astype(float);noise=variance[ti,columns].astype(float)
        total=noise[:,None]+scales[None,:]**2
        logp=-.5*(np.log(total)+value[:,None]**2/total)+np.log(np.maximum(mix,1e-12))
        responsibility=np.exp(logp-logsumexp(logp,axis=1)[:,None])
        gain=(responsibility*scales[None,:]**2/total).sum(axis=1)
        out[ti,columns]=value*gain;priors.append(mix)
    return out,np.asarray(priors)
def predict_modes(data,train,query_context):
    cube=collapse(data,train);outputs={};details={}
    for kind in ['uniform','rank']:
        weights,rho=context_weights(data,query_context,train,kind=='rank')
        mean,var,valid=combine(cube,weights)
        outputs[kind]=mean;details[kind]={'weights':weights.tolist(),'rho':rho.tolist()}
        if kind=='rank':
            outputs['adaptive'],mixture=adaptive_shrink(mean,var,valid)
            details['adaptive']={'mixture':mixture.tolist(),**details[kind]}
            m,v,ok=combine(cube,weights,True)
            outputs['adaptive_kd'],mix=adaptive_shrink(m,v,ok)
            details['adaptive_kd']={'mixture':mix.tolist(),**details[kind]}
    return outputs,details
def bh(p,valid):
    q=np.ones_like(p)
    cols=np.flatnonzero(valid)
    if len(cols):
        order=np.argsort(p[cols],kind='stable');rank=np.arange(1,len(cols)+1)
        result=np.minimum.accumulate((p[cols][order]*len(cols)/rank)[::-1])[::-1]
        q[cols[order]]=np.minimum(result,1.)
    return q
def summary_score(pred,truth,error,q,valid):
    real=valid&(q<.05);p=2*norm.sf(np.abs(pred)/np.maximum(error,.05))
    called=valid&(bh(p,valid)<.05);correct=called&(np.sign(pred)==np.sign(truth))
    budget=max(int(real.sum()),int(called.sum()));union=real|called
    mse=float(((pred-truth)[valid]**2).mean());zero=float((truth[valid]**2).mean())
    fidelity=float(correct.sum()/budget) if budget else 0.
    jac=float((real&called).sum()/union.sum()) if union.any() else 1.
    ratio=mse/max(zero,1e-8)
    return {'mse':mse,'zero_mse':zero,'mse_ratio':ratio,'fidelity_proxy':fidelity,'jaccard_proxy':jac,
        'predicted_calls':int(called.sum()),'reference_calls':int(real.sum()),
        'sign_reference_DE':float((np.sign(pred[real])==np.sign(truth[real])).mean()) if real.any() else np.nan,
        'selection_objective':.5*fidelity+.25*jac+.25*(1-min(ratio,2.))}

