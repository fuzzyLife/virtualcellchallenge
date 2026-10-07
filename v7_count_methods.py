"""V7 integer emission: original V4, reduced intervention noise, or partial moments."""
import numpy as np
from scipy import sparse
import build_submission as build
from v6_methods import moments
from build_v5_submission import project_bulk
def round_unbiased(value,rng):
    return np.floor(value+rng.random(value.shape)).astype(np.int64)
def emit(x,rows,effect,own,seed,generator='v4',atlas=None,prepared=None):
    if prepared is None:
        depth,p,b,_=moments(x)
        prepared=(depth,p,b,np.asarray(x.mean(axis=0)).ravel())
    depth,p,b,mean_raw=prepared
    factor=np.exp2(np.clip(effect,-3,3))
    if own is not None:factor[own]=.1
    if generator=='v4':
        return build.simulate_counts_zero_aware(x,rows,factor,np.random.default_rng(seed),mean_raw,depth.mean())
    if generator=='cal10':
        initial=build.simulate_counts_zero_aware(x,rows,factor,np.random.default_rng(seed),mean_raw,depth.mean())
        libs,pi,bi,_=moments(initial)
        wanted_p=p*factor;wanted_p/=wanted_p.sum()
        wanted_b=b*factor;wanted_b/=wanted_b.sum()
        desired_p=.9*pi+.1*wanted_p;desired_b=.9*bi+.1*wanted_b
        projected=project_bulk(desired_p,desired_b,libs)
        template=initial.toarray().astype(float)/libs[:,None]
        count=atlas.dual_moment_counts(template,desired_p,desired_b,depths=libs.astype(np.int64),seed=int(seed),iterations=200)
        np.testing.assert_array_equal(count.sum(axis=1),libs)
        _,pf,bf,_=moments(count)
        if np.abs(pf-desired_p).sum()>.03 or np.abs(bf-projected).sum()>.002:raise ValueError('Moment calibration failed')
        return sparse.csr_matrix(count)
    if generator!='low_noise':raise ValueError('Unknown count generator')
    rng=np.random.default_rng(seed);mix=np.sqrt(.5)
    original=sparse.csr_matrix(x)[rows].astype(np.int64)
    result=original.copy();idx=result.indices
    lower=factor[idx]<1
    expectation=result.data[lower]*factor[idx[lower]]
    draw=rng.binomial(result.data[lower],factor[idx[lower]])
    result.data[lower]=round_unbiased((1-mix)*expectation+mix*draw,rng)
    library=np.asarray(original.sum(axis=1)).ravel()/depth.mean()
    increased=np.flatnonzero(factor>1)
    for start in range(0,len(increased),128):
        cols=increased[start:start+128]
        rates=library[:,None]*(mean_raw[cols]*(factor[cols]-1))[None,:]
        draw=rng.poisson(rates);add=round_unbiased((1-mix)*rates+mix*draw,rng)
        rr,cc=np.nonzero(add)
        if len(rr):result+=sparse.csr_matrix((add[rr,cc],(rr,cols[cc])),shape=result.shape)
    result.eliminate_zeros();result.sort_indices()
    if (result.data<0).any() or result.data.max(initial=0)>np.iinfo(np.int32).max:raise ValueError('Invalid count emission')
    return result.astype(np.int32)

