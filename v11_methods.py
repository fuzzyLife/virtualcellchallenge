"""V11 consensus-supported response coverage with isolated population emission dials."""
from pathlib import Path
import numpy as np
from scipy import sparse
from v7_response_methods import collapse,combine,adaptive_shrink
from v8_methods import old_response,gated
from v6_methods import moments
import build_submission as build
ROOT=build.ROOT;OUT=ROOT/'experiments/v11_consensus_coverage';A=OUT/'analysis'
SEEDS=[20281008,20281031];SEED=20261011
SETTINGS={
 'v10':{'response':'dose','amplitude':1.,'noise':.5,'bulk_amplitude':1.,'eligible':True},
 'v4_compact':{'response':'v4','amplitude':1.,'noise':.5,'bulk_amplitude':1.,'eligible':False},
 'dose125':{'response':'dose','amplitude':1.25,'noise':.5,'bulk_amplitude':1.,'eligible':True},
 'consensus20':{'response':'consensus20','amplitude':1.,'noise':.5,'bulk_amplitude':1.,'eligible':True},
 'consensus40':{'response':'consensus40','amplitude':1.,'noise':.5,'bulk_amplitude':1.,'eligible':True},
 'noise35':{'response':'dose','amplitude':1.,'noise':.35,'bulk_amplitude':1.,'eligible':True},
 'noise65':{'response':'dose','amplitude':1.,'noise':.65,'bulk_amplitude':1.,'eligible':True},
 'consensus20_noise65':{'response':'consensus20','amplitude':1.,'noise':.65,'bulk_amplitude':1.,'eligible':True},
 'bulk085':{'response':'dose','amplitude':1.,'noise':.5,'bulk_amplitude':.85,'eligible':True},
}


def consensus_gate(y,var,valid,kd,weights,mean,variance):
    confidence=np.clip(.04/np.maximum(var,.0025),.1,4.)
    w=np.asarray(weights)[:,None,None]*confidence*valid;den=w.sum(axis=0)
    signed=(w*np.sign(y)).sum(axis=0)
    agreement=np.divide(abs(signed),den,out=np.zeros_like(mean,dtype=float),where=den>0)
    # Rescue only directions supported by at least two independently measured
    # families. Contradictory sources and high uncertainty keep the V10 estimate.
    support=valid.sum(axis=0)
    gate=np.clip((agreement-.5)/.5,0,1)*np.minimum(abs(mean)/np.sqrt(variance),1.)*(support>=2)
    return gate,agreement,support


def response_bank(data,excluded=None):
    train=np.ones(len(data['truth']),bool) if excluded is None else data['row_families']!=excluded
    y,var,valid,kd=collapse(data,train);weights=np.array([2.,.25,2.,1.,1.])
    gain=np.clip(.9/np.maximum(kd,.35),.8,1.8);y*=gain[:,:,None];var*=gain[:,:,None]**2
    mean,variance,ok=combine((y,var,valid,kd),weights)
    eb,prior=adaptive_shrink(mean,variance,ok)
    gate,agreement,support=consensus_gate(y,var,valid,kd,weights,mean,variance)
    bank={'dose':np.clip(eb,-3,3).astype(np.float32),'v4':old_response(line=excluded)}
    for alpha in [.2,.4]:
        bank[f'consensus{int(alpha*100)}']=np.clip(eb+alpha*gate*(mean-eb),-3,3).astype(np.float32)
    return bank,{'excluded_family':excluded,'training_families':sorted(set(data['row_families'][train].tolist())),
        'train_rows':np.flatnonzero(train).tolist(),'family_weights':weights.tolist(),
        'dose_gain_bounds':[.8,1.8],'consensus_min_families':2,'consensus_agreement_floor':.5,
        'rescued_fraction':float((gate>0).mean()),'mean_rescue_gate':float(gate.mean())}


def response_for(bank,setting,query,own,target_index):
    return gated(bank[setting['response']][target_index]*setting['amplitude'],query,own).astype(np.float32)


def emit(x,rows,effect,own,seed,setting,prepared=None,atlas=None):
    if prepared is None:depth,p,b,_=moments(x)
    else:depth,p,b=prepared[:3]
    rows=np.asarray(rows,dtype=int);libs=np.asarray(depth)[rows].astype(np.int64)
    # Indexing must describe x, even when x is a sampled validation template.
    if len(depth)!=x.shape[0]:raise ValueError('Depth axis must match the supplied cell matrix')
    np.testing.assert_array_equal(libs,np.asarray(x[rows].sum(axis=1)).ravel())
    factor=np.exp2(np.clip(effect,-3,3));bulk_factor=np.exp2(np.clip(effect*setting['bulk_amplitude'],-3,3))
    if own is not None:factor[own]=bulk_factor[own]=.1
    wanted_p=p*factor;wanted_p/=wanted_p.sum();wanted_b=b*bulk_factor;wanted_b/=wanted_b.sum()
    noise=float(setting['noise']);assert 0<=noise<=1
    rng=np.random.default_rng(seed);expected=libs[:,None]*wanted_p[None,:]
    noisy=rng.poisson(expected);template=(1-noise)*expected+noise*noisy
    template/=template.sum(axis=1,keepdims=True)
    if atlas is None:
        from evaluate_atlas_transfer import load_atlas_model
        atlas,_=load_atlas_model()
    count=atlas.dual_moment_counts(template,wanted_p,wanted_b,depths=libs,seed=int(seed),iterations=200)
    np.testing.assert_array_equal(count.sum(axis=1,dtype=np.int64),libs)
    if not np.issubdtype(count.dtype,np.integer) or (count<0).any():raise ValueError('Invalid raw counts')
    realized_p=np.mean(count/libs[:,None],axis=0);realized_b=count.sum(axis=0,dtype=np.int64)/libs.sum()
    diag={'noise_amplitude':noise,'bulk_response_amplitude':setting['bulk_amplitude'],
          'mean_probability_L1':float(abs(realized_p-wanted_p).sum()),'bulk_probability_L1':float(abs(realized_b-wanted_b).sum())}
    if diag['mean_probability_L1']>.03 or diag['bulk_probability_L1']>.031:raise ValueError('Population transport constraints failed')
    return sparse.csr_matrix(count),diag
