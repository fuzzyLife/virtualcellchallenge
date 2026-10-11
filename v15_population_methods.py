"""V15 source-treated marginal changes applied to ranked query-control cells."""
import json,gc
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.stats import rankdata
from v15_ensemble_methods import ROOT,A,FAMILIES,cube
from fit_v7_response import raw_inputs
from v6_methods import moments
from v11_methods import emit as compact,SETTINGS
import build_submission as build
Q=np.array([0.,.25,.5,.75,1.])

def source_quantiles(data):
    folder=A/'distribution_sources';folder.mkdir(parents=True,exist_ok=True);records=[]
    for fi,family in enumerate(FAMILIES):
        if family not in ['H1','HCT116','HEK293T']:continue
        path=folder/(str(family)+'.npz');side=path.with_suffix('.json')
        inputs={'source_cache':build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),'code':build.sha256(Path(__file__))}
        if path.exists() and side.exists():
            rec=json.loads(side.read_text());assert rec['inputs']==inputs and build.sha256(path)==rec['sha256'];records.append(rec);continue
        x,meta,measured,provenance=raw_inputs(str(family),data['genes']);controls=np.flatnonzero(meta.gene_target.eq('Non-Targeting'))
        depths,cp,_,_=moments(x[controls]);split=float(np.median(depths));qctrl=[];counts=[]
        for bi in range(2):
            available=controls[(depths<=split) if bi==0 else (depths>split)]
            ids=np.random.default_rng(2026101501+fi*10+bi).choice(available,min(len(available),800),replace=False)
            libs=np.asarray(x[ids].sum(1)).ravel();relative=x[ids].toarray().astype(float)/libs[:,None]
            relative=np.divide(relative,cp[None,:],out=np.zeros_like(relative),where=cp[None,:]>0)
            qctrl.append(np.quantile(relative,Q,axis=0).T.astype(np.float32));counts.append(len(ids))
        delta=np.zeros((300,2,len(data['genes']),5),np.float32);valid=np.zeros((300,len(data['genes'])),bool);n=np.zeros(300,int)
        for ti,target in enumerate(data['targets']):
            ids=np.flatnonzero(meta.gene_target.eq(target));n[ti]=len(ids)
            if len(ids)<80:continue
            treated=x[ids];libs,tp,_,_=moments(treated);mask=measured&(cp*1e6>5)&(tp*1e6>5)
            valid[ti]=mask
            for bi in range(2):
                use=(libs<=split) if bi==0 else (libs>split)
                if use.sum()<25:valid[ti]=False;continue
                relative=treated[use].toarray().astype(float)/libs[use,None]
                relative=np.divide(relative,tp[None,:],out=np.zeros_like(relative),where=tp[None,:]>0)
                qt=np.quantile(relative,Q,axis=0).T
                # Quantile changes are relative to each population's own mean; response means are separate.
                delta[ti,bi]=np.clip(qt-qctrl[bi],-4,4)
            delta[ti,:,~valid[ti],:]=0
            if (ti+1)%75==0:print('V15 distribution source',family,ti+1,flush=True)
        np.savez_compressed(path,delta=delta,valid=valid,n=n,control_quantiles=np.array(qctrl),depth_split=split)
        rec={'family':str(family),'inputs':inputs,'sha256':build.sha256(path),'eligible_targets':int(valid.any(1).sum()),'control_quantile_cell_counts':counts,'source':provenance,
             'policy':'Raw source controls/treatments only; normalized within-population means, two source-control depth bins, five quantiles, minimum 80 treatments and 25 per depth bin; source library totals projected to competition gene axis.'}
        build.write_json(side,rec);records.append(rec);del x,delta;gc.collect();print('V15 distribution source complete',family,flush=True)
    return records

def distribution_model(data,excluded):
    y,var,measured,own=cube(data,excluded);result=np.zeros((300,2,len(data['genes']),5),np.float32);den=np.zeros((300,len(data['genes'])),float);support=np.zeros_like(den,int)
    for fi,f in enumerate(FAMILIES):
        if f in excluded or f not in ['H1','HCT116','HEK293T']:continue
        path=A/'distribution_sources'/(str(f)+'.npz');side=json.loads(path.with_suffix('.json').read_text());assert side['sha256']==build.sha256(path)
        with np.load(path,allow_pickle=False) as z:
            valid=z['valid']&measured[fi];w=np.array([2.,.25,2.,1.,1.])[fi]*np.clip(.04/np.maximum(var[fi],.0025),.1,4.)*valid*np.minimum(z['n'][:,None]/200,1.)
            result+=z['delta']*w[:,None,:,None];den+=w;support+=valid
    result=np.divide(result,den[:,None,:,None],out=np.zeros_like(result),where=den[:,None,:,None]>0)
    result[np.arange(300),:,own,:]=0
    return result,den>0,support

def transport_initial(template,libs,p,wanted,delta,mask,seed):
    relative=template.toarray().astype(float)/libs[:,None]
    relative=np.divide(relative,p[None,:],out=np.zeros_like(relative),where=p[None,:]>0)
    rng=np.random.default_rng(seed);bins=libs>np.median(libs)
    # Randomized ordering of exact ties distributes zero cells across their empirical CDF interval.
    for bi in range(2):
        ids=np.flatnonzero(bins==bool(bi))
        if not len(ids):continue
        r=relative[ids];tie=rng.random(r.shape)
        order=np.lexsort((tie,r),axis=0);ranks=np.empty_like(order);np.put_along_axis(ranks,order,np.arange(len(ids))[:,None],axis=0)
        u=(ranks+.5)/len(ids);position=np.minimum((u*4).astype(int),3);frac=u*4-position
        left=np.take_along_axis(delta[bi][None,:,:],position[:,:,None],axis=2)[...,0]
        right=np.take_along_axis(delta[bi][None,:,:],(position+1)[:,:,None],axis=2)[...,0]
        change=left*(1-frac)+right*frac
        relative[ids]=np.where(mask[None,:],np.maximum(r+.5*change,0),r)
    return libs[:,None]*wanted[None,:]*relative

def emit(x,rows,effect,own,seed,kind,delta=None,mask=None,prepared=None,atlas=None):
    if kind=='compact':return compact(x,rows,effect,own,seed,SETTINGS['consensus20_noise65'],prepared=prepared,atlas=atlas)
    depth,p,b=prepared[:3] if prepared is not None else moments(x)[:3]
    rows=np.asarray(rows,int);libs=np.asarray(depth)[rows].astype(np.int64)
    assert len(depth)==x.shape[0];np.testing.assert_array_equal(libs,np.asarray(x[rows].sum(1)).ravel())
    factor=np.exp2(np.clip(effect,-3,3))
    if own is not None:factor[own]=.1
    wp=p*factor;wp/=wp.sum();wb=b*factor;wb/=wb.sum()
    linked=transport_initial(x[rows],libs,p,wp,delta,mask,seed+415)
    rng=np.random.default_rng(seed);expected=libs[:,None]*wp[None,:];initial=.35*expected+.65*rng.poisson(expected)
    alpha=.1 if kind.endswith('10') else .25
    # Add empirical state/marginal variation only on measured, source-supported output genes.
    initial[:,mask]=(1-alpha)*initial[:,mask]+alpha*linked[:,mask]
    initial/=initial.sum(1,keepdims=True)
    if atlas is None:
        from evaluate_atlas_transfer import load_atlas_model
        atlas,_=load_atlas_model()
    counts=atlas.dual_moment_counts(initial,wp,wb,depths=libs,seed=int(seed),iterations=200)
    np.testing.assert_array_equal(counts.sum(1,dtype=np.int64),libs);assert np.issubdtype(counts.dtype,np.integer) and (counts>=0).all()
    rp=np.mean(counts/libs[:,None],axis=0);rb=counts.sum(0,dtype=np.int64)/libs.sum()
    diag={'mean_probability_L1':float(abs(rp-wp).sum()),'bulk_probability_L1':float(abs(rb-wb).sum()),'distribution_blend':alpha,'distribution_supported_genes':int(mask.sum())}
    if diag['mean_probability_L1']>.03 or diag['bulk_probability_L1']>.031:raise ValueError('Transport moment constraints failed')
    return sparse.csr_matrix(counts),diag
