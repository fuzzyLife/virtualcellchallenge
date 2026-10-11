"""V15: nested held-family ensemble of pooled, robust and retrieval responses."""
from pathlib import Path
import json,gc
import numpy as np
from v7_response_methods import FAMILIES,collapse,adaptive_shrink,summary_score
from v8_methods import gated
from analyze_target_rank_correlations import pairwise_spearman
import build_submission as build
ROOT=build.ROOT;OUT=ROOT/'experiments/v15_weighted_transport';A=OUT/'analysis'
EXPERTS=['v13','student','specialist']
GRID=[(1.,0.,0.),(.75,.25,0.),(.5,.5,0.),(.75,0.,.25),(.5,.25,.25),(.5,0.,.5)]

def cube(data,excluded):
    own=np.array([np.flatnonzero(data['genes']==t)[0] for t in data['targets']])
    d=dict(data);d['knockdown']=data['knockdown'].copy()
    d['truth']=np.where(data['valid'],data['truth'],0);d['errors']=np.where(data['valid'],data['errors'],1e3)
    p=data['profiles'][data['row_contexts'],own[data['row_targets']]]
    d['knockdown'][~(np.isfinite(p)&(p>5))]=.9
    train=~np.isin(data['row_families'],excluded)
    y,var,valid,kd=collapse(d,train);gain=np.clip(.9/np.maximum(kd,.35),.8,1.8)
    y*=gain[:,:,None];var*=gain[:,:,None]**2
    return y,var,valid,own

def profiles(data):
    out=[]
    for f in FAMILIES:
        ids=np.unique(data['row_contexts'][data['row_families']==f]);p=data['profiles'][ids];m=np.isfinite(p);n=m.sum(0)
        out.append(np.divide(np.where(m,p,0).sum(0),n,out=np.full(len(data['genes']),np.nan),where=n>0))
    return np.asarray(out)

def pool(y,var,valid,weights,own,robust=False):
    y=np.where(valid,y,0);var=np.where(valid,var,1e6)
    initial=weights[:,None,None]*np.clip(.04/np.maximum(var,.0025),.1,4.)*valid
    w=initial
    if robust:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            location=np.nan_to_num(np.nanmedian(np.where(valid,y,np.nan),axis=0))
            noise=np.nan_to_num(np.nanmedian(np.where(valid,var,np.nan),axis=0),nan=1.)
        for _ in range(5):
            w=initial*4/(3+(y-location[None,:])**2/(var+np.maximum(.25**2,noise)[None,:]))
            den=w.sum(0);wn=np.divide(w,den[None,:],out=np.zeros_like(w),where=den[None,:]>0);location=(wn*y).sum(0)
    den=w.sum(0);wn=np.divide(w,den[None,:],out=np.zeros_like(w),where=den[None,:]>0);ok=den>0
    mean=(wn*y).sum(0);variance=np.maximum((wn*wn*var).sum(0)+.25*(wn*(y-mean[None,:])**2).sum(0),.0025)
    mean[~ok]=0;variance[~ok]=1e6;eb,_=adaptive_shrink(mean,variance,ok)
    agreement=np.divide(abs((w*np.sign(y)).sum(0)),den,out=np.zeros_like(mean),where=den>0)
    gate=np.clip((agreement-.5)/.5,0,1)*np.minimum(abs(mean)/np.sqrt(variance),1.)*(valid.sum(0)>=2)
    effect=np.clip(eb+.2*gate*(mean-eb),-3,3).astype(np.float32);effect[np.arange(len(own)),own]=0
    return effect

def cache_experts(data,excluded):
    folder=A/'response_cache';folder.mkdir(parents=True,exist_ok=True)
    key='_'.join(sorted(excluded)) or 'all';path=folder/(key+'.npz');side=path.with_suffix('.json')
    inputs={'source':build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),'code':build.sha256(Path(__file__)),'excluded':sorted(excluded)}
    if path.exists() and side.exists():
        record=json.loads(side.read_text());assert record['inputs']==inputs and build.sha256(path)==record['sha256']
        with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
    y,var,valid,own=cube(data,excluded);weights=np.array([2.,.25,2.,1.,1.])
    arrays={'v13':pool(y,var,valid,weights,own),'student':pool(y,var,valid,weights,own,True),'own':own}
    np.savez_compressed(path,**arrays);build.write_json(side,{'inputs':inputs,'sha256':build.sha256(path),'training_families':[str(f) for f in FAMILIES if f not in excluded]})
    print('V15 cached response experts:',key,flush=True);return arrays

def source_reliability(data,excluded):
    y,var,valid,own=cube(data,excluded);axis=~np.isin(data['genes'],data['targets']);loss=np.full((5,5),np.nan)
    for i in range(5):
        if FAMILIES[i] in excluded:continue
        for j in range(5):
            if i==j or FAMILIES[j] in excluded:continue
            m=valid[i]&valid[j]&axis[None,:]
            # Equal target votes; normalized transfer loss is bounded and avoids near-zero denominators.
            numerator=np.where(m,(np.clip(y[i],-3,3)-np.clip(y[j],-3,3))**2,0).sum(1)
            denominator=np.where(m,y[j]**2+var[j],0).sum(1)
            ok=(m.sum(1)>=100)&(denominator>0)
            if ok.any():loss[i,j]=np.median(np.clip(numerator[ok]/denominator[ok],0,4))
    reliability=np.ones(5)
    for i in range(5):
        v=loss[i,np.isfinite(loss[i])]
        if len(v):reliability[i]=float(np.clip(np.exp(1-np.mean(v)),.25,4.))
    reliability[np.isin(FAMILIES,excluded)]=0
    return reliability,loss

def specialist(data,excluded,query,reliability):
    y,var,valid,own=cube(data,excluded);p=profiles(data);axis=~np.isin(data['genes'],data['targets']);rho=np.full(5,np.nan)
    for i in range(5):
        if FAMILIES[i] in excluded:continue
        r,n,status=pairwise_spearman(query[axis],p[i,axis],min_pairs=100)
        if status=='ok':rho[i]=r
    weights=reliability.copy();ok=np.isfinite(rho)
    if ok.any():weights[ok]*=np.exp(2*(rho[ok]-np.max(rho[ok])))
    weights[~ok]=0
    return pool(y,var,valid,weights,own),{'weights':weights.tolist(),'rho':rho.tolist()}

def gate_bank(bank,query,own):
    p=np.nan_to_num(query,nan=0.)
    return {k:np.asarray([gated(v.copy(),p,int(o)) for v,o in zip(a,own)],dtype=np.float32) for k,a in bank.items()}

def fit_ensemble(data,excluded=()):
    excluded=tuple(sorted(excluded));folds=[str(f) for f in FAMILIES if f not in excluded];score=np.zeros((len(folds),len(GRID)))
    for fi,f in enumerate(folds):
        inner=tuple(sorted((*excluded,f)));base=cache_experts(data,inner);reliability,_=source_reliability(data,inner)
        rows=np.flatnonzero((data['row_families']==f)&~np.isin(data['targets'][data['row_targets']],json.loads((A/'plan.json').read_text())['reserved_targets']))
        context_scores=[]
        for ci in np.unique(data['row_contexts'][rows]):
            sp,_=specialist(data,inner,data['profiles'][ci],reliability)
            bank=gate_bank({'v13':base['v13'],'student':base['student'],'specialist':sp},data['profiles'][ci],base['own']);values=[]
            for ri in rows[data['row_contexts'][rows]==ci]:
                ti=data['row_targets'][ri];mask=data['valid'][ri].copy();mask[base['own'][ti]]=False
                pred=np.asarray([bank[n][ti] for n in EXPERTS]);one=[]
                for weight in GRID:
                    effect=np.einsum('i,ig->g',weight,pred);one.append(summary_score(effect,data['truth'][ri],data['errors'][ri],data['q'][ri],mask)['selection_objective'])
                values.append(one)
            if values:context_scores.append(np.mean(values,axis=0))
        score[fi]=np.mean(context_scores,axis=0) if context_scores else np.nan
        print('V15 inner ensemble calibration:',excluded,f,flush=True)
    means=np.nanmean(score,axis=0);best=int(np.nanargmax(means));base=cache_experts(data,excluded);reliability,loss=source_reliability(data,excluded)
    return base,{'excluded':list(excluded),'weights':list(GRID[best]),'inner_families':folds,'inner_objectives':score.tolist(),'inner_mean_objectives':means.tolist(),'source_reliability':reliability.tolist(),'source_pair_losses':np.nan_to_num(loss,nan=-1).tolist()}

def predict(data,base,model,query):
    sp,diag=specialist(data,model['excluded'],query,np.asarray(model['source_reliability']))
    bank=gate_bank({'v13':base['v13'],'student':base['student'],'specialist':sp},query,base['own'])
    bank['student50']=.5*bank['v13']+.5*bank['student']
    bank['ensemble']=sum(w*bank[n] for w,n in zip(model['weights'],EXPERTS)).astype(np.float32)
    return bank,diag
