"""V14: masked, within-target ridge interactions on training-only control ranks."""
from pathlib import Path
import numpy as np
from scipy.stats import rankdata
from v7_response_methods import FAMILIES, collapse
from v13_quality_methods import banks
from v8_methods import gated
import build_submission as build
ROOT=build.ROOT
A=ROOT/'experiments/v14_context_interactions/analysis'
NAMES=['v13','context2','context8','shuffled2','amplitude110']


def family_profiles(data):
    result=[]
    for family in FAMILIES:
        ids=np.unique(data['row_contexts'][data['row_families']==family])
        if not len(ids):result.append(np.full(len(data['genes']),np.nan));continue
        p=data['profiles'][ids];valid=np.isfinite(p);den=valid.sum(0)
        result.append(np.divide(np.where(valid,p,0).sum(0),den,out=np.full(p.shape[1],np.nan),where=den>0))
    return np.asarray(result)


def feature_fit(profiles,training,genes,targets):
    axis=np.flatnonzero(np.isfinite(profiles[training]).all(0)&~np.isin(genes,targets))
    if len(axis)<100:raise ValueError('Too few shared non-panel control genes')
    ranked=np.array([rankdata(p[axis],method='average')/len(axis) for p in profiles[training]])
    center=ranked.mean(0);_,_,vt=np.linalg.svd(ranked-center,full_matrices=False)
    basis=vt[:2];scores=(ranked-center)@basis.T
    scale=np.maximum(np.sqrt(np.mean(scores**2,axis=0)),1e-6)
    scores/=scale
    radius=max(float(np.linalg.norm(scores,axis=1).max())*1.5,1.)
    return {'axis':axis,'center':center,'basis':basis,'scale':scale,'radius':radius},scores


def feature_query(model,profile):
    values=np.asarray(profile)[model['axis']];valid=np.isfinite(values)
    coverage=float(valid.mean())
    if coverage<.7:return np.zeros(2),{'coverage':coverage,'abstain':True,'clipped':False}
    rank=model['center'].copy();rank[valid]=rankdata(values[valid],method='average')/valid.sum()
    z=((rank-model['center'])@model['basis'].T)/model['scale']
    length=float(np.linalg.norm(z));clipped=length>model['radius']
    if clipped:z*=model['radius']/length
    return z,{'coverage':coverage,'abstain':False,'clipped':bool(clipped)}


def ridge_interactions(y,var,valid,z,penalty):
    """Fit centered family effects; invalid entries abstain, never vote zero.

    y: families x targets x genes. Every family has one vote before uncertainty
    weighting, including CD4's collapsed conditions. Return target/gene slopes.
    """
    w=np.where(valid,np.clip(.04/np.maximum(var,.0025),.1,4.),0.)
    den=w.sum(0);support=valid.sum(0)
    normalized=np.divide(w,den[None,:],out=np.zeros_like(w,dtype=float),where=den[None,:]>0)
    center=np.einsum('ftg,fk->tgk',normalized,z)
    average=np.sum(normalized*np.where(valid,y,0),axis=0)
    xt=z[:,None,None,:]-center[None,:]
    # Normalize total confidence so shrinkage is not a hidden cell-count gain.
    weights=normalized*support[None,:]
    gram=np.einsum('ftg,ftgk,ftgl->tgkl',weights,xt,xt)
    gram+=penalty*np.eye(z.shape[1])[None,None,:]
    rhs=np.einsum('ftg,ftgk,ftg->tgk',weights,xt,np.where(valid,y-average[None,:],0))
    beta=np.linalg.solve(gram,rhs[...,None])[...,0]
    beta[support<3]=0
    assert np.isfinite(beta).all()
    return beta.astype(np.float32),center.astype(np.float32),support


def fit(data,excluded):
    if excluded not in FAMILIES:raise ValueError('Unknown held-out family')
    training=np.flatnonzero(FAMILIES!=excluded)
    train=data['row_families']!=excluded
    adjusted=dict(data);adjusted['knockdown']=data['knockdown'].copy()
    own=np.array([np.flatnonzero(data['genes']==t)[0] for t in data['targets']])
    control=data['profiles'][data['row_contexts'],own[data['row_targets']]]
    known=np.isfinite(control)&(control>5);adjusted['knockdown'][~known]=.9
    y,var,valid,kd=collapse(adjusted,train)
    gain=np.clip(.9/np.maximum(kd,.35),.8,1.8);y*=gain[:,:,None];var*=gain[:,:,None]**2
    baseline,record=banks(data,excluded)
    base=baseline['unknown_kd']['consensus20'].copy()
    profiles=family_profiles(data)
    features,z=feature_fit(profiles,training,data['genes'],data['targets'])
    result={'base':base,'features':features,'own':own,'excluded':excluded,'training_families':FAMILIES[training].tolist(),'slopes':{},'centers':{}}
    for name,penalty,shuffled in [('context2',2.,False),('context8',8.,False),('shuffled2',2.,True)]:
        design=z.copy()
        if shuffled:design=design[np.random.default_rng(2026101403).permutation(len(design))]
        beta,center,support=ridge_interactions(y[training],var[training],valid[training],design,penalty)
        beta[np.arange(len(own)),own]=0
        result['slopes'][name]=beta;result['centers'][name]=center
    result['training_rows']=np.flatnonzero(train).tolist()
    result['feature_genes']=len(features['axis'])
    assert excluded not in result['training_families']
    return result


def predict(model,profile,apply_gate=True):
    z,diagnostic=feature_query(model['features'],profile)
    result={'v13':model['base'].copy(),'amplitude110':np.clip(model['base']*1.1,-3,3)}
    for name,beta in model['slopes'].items():
        correction=np.sum(beta*(z-model['centers'][name]),axis=-1)
        if diagnostic['abstain']:correction[:]=0
        result[name]=np.clip(model['base']+.5*np.clip(correction,-.5,.5),-3,3).astype(np.float32)
    for name,value in result.items():
        value[np.arange(len(model['own'])),model['own']]=0
        if apply_gate:
            result[name]=np.array([gated(row.copy(),np.nan_to_num(profile,nan=0.),int(own)) for row,own in zip(value,model['own'])],dtype=np.float32)
    return result,diagnostic
