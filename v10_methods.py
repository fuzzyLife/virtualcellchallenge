"""V10 measured-intervention pooling and controlled population count transport.

No challenge treated labels enter either response fitting or emission.
"""
from pathlib import Path
import numpy as np
from scipy import sparse
from v7_response_methods import collapse, combine, adaptive_shrink
from v8_methods import old_response
from v6_methods import moments
import build_submission as build
ROOT=build.ROOT
OUT=ROOT/'experiments/v10_population_transport'
A=OUT/'analysis'
SEEDS=[20271008,20271027]
SETTINGS={
 'v4_legacy':('v4','legacy'),
 'v4_compact':('v4','compact'),
 'pooled_compact':('pooled','compact'),
 'dose_compact':('dose','compact'),
 'invariant_compact':('invariant','compact'),
 'dose_selective':('dose','selective'),
}

def response_bank(data,excluded=None,context=None):
    train=np.ones(len(data['truth']),bool) if excluded is None else data['row_families']!=excluded
    y,var,valid,kd=collapse(data,train)
    weights=np.array([2.,.25,2.,1.,1.])
    bank={'v4':old_response(line=excluded,context=context)}
    for mode in ['pooled','dose','invariant']:
        yy=y.copy();vv=var.copy()
        if mode!='pooled':
            gain=np.clip(.9/np.maximum(kd,.35),.8,1.8)
            yy*=gain[:,:,None];vv*=gain[:,:,None]**2
        if mode=='invariant':
            # Subtract part of the source-wide shift: perturbation-specific signal
            # is less sensitive to source batch/platform and shared stress axes.
            den=valid.sum(axis=1)
            common=np.divide((yy*valid).sum(axis=1),den,out=np.zeros_like(den,dtype=float),where=den>0)
            yy-=.5*common[:,None,:]
        mean,variance,ok=combine((yy,vv,valid,kd),weights)
        value,prior=adaptive_shrink(mean,variance,ok)
        bank[mode]=np.clip(value,-3,3).astype(np.float32)
    return bank,{'excluded_family':excluded,'training_families':sorted(set(data['row_families'][train].tolist())),
                 'family_weights':weights.tolist(),'dose_gain_bounds':[.8,1.8],
                 'source_center_fraction':.5,'train_rows':np.flatnonzero(train).tolist()}

def emit(x,rows,effect,own,seed,generator='compact',prepared=None,atlas=None):
    if prepared is None:
        depth,p,b,_=moments(x);raw=np.asarray(x.mean(axis=0)).ravel()
    else:depth,p,b,raw=prepared
    factor=np.exp2(np.clip(effect,-3,3))
    if own is not None:factor[own]=.1
    if generator=='legacy':
        return build.simulate_counts_zero_aware(x,rows,factor,np.random.default_rng(seed),raw,depth.mean()),{'generator':generator}
    if generator not in ['compact','selective']:raise ValueError(generator)
    if atlas is None:
        from evaluate_atlas_transfer import load_atlas_model
        atlas,_=load_atlas_model()
    libs=depth[rows].astype(np.int64)
    wanted_p=p*factor;wanted_p/=wanted_p.sum()
    wanted_b=b*factor;wanted_b/=wanted_b.sum()
    rng=np.random.default_rng(seed)
    # Independent synthetic profiles with half-Poisson amplitude; use the
    # real observed depth distribution, but not each cell's expression vector.
    expected=libs[:,None]*wanted_p[None,:]
    noisy=rng.poisson(expected)
    template=.5*expected+.5*noisy
    template/=template.sum(axis=1,keepdims=True)
    if generator=='selective':
        original=sparse.csr_matrix(x)[rows].toarray().astype(float)*factor
        original/=original.sum(axis=1,keepdims=True)
        precision=np.clip(np.abs(effect)/.25,0,1)
        if own is not None:precision[own]=1
        template=precision[None,:]*template+(1-precision[None,:])*original
        template/=template.sum(axis=1,keepdims=True)
    count=atlas.dual_moment_counts(template,wanted_p,wanted_b,depths=libs,seed=int(seed),iterations=200)
    np.testing.assert_array_equal(count.sum(axis=1,dtype=np.int64),libs)
    realized_p=np.mean(count/libs[:,None],axis=0)
    realized_b=count.sum(axis=0,dtype=np.int64)/libs.sum()
    diag={'generator':generator,'mean_probability_L1':float(np.abs(realized_p-wanted_p).sum()),
          'bulk_probability_L1':float(np.abs(realized_b-wanted_b).sum()),
          'dispersion_amplitude':.5,'response_selected_fraction':float((np.abs(effect)>=.25).mean())}
    if diag['mean_probability_L1']>.03 or diag['bulk_probability_L1']>.031:raise ValueError('Population transport constraints failed')
    return sparse.csr_matrix(count),diag
