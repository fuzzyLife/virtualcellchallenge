"""V68: depth-balanced sampling of real control cells to model a response."""
from pathlib import Path
import argparse
import numpy as np
from scipy import sparse
from v57_pie import read,write,sha

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v68_population_balancing';A=OUT/'analysis'
RIDGE=1e-4


def strata(libraries,draw,bins=8):
    n=len(libraries);assert n>=draw
    order=np.argsort(libraries,kind='stable');group=np.empty(n,int)
    group[order]=np.arange(n)*min(bins,n)//n
    sizes=np.bincount(group);desired=sizes*draw/n;quota=np.floor(desired).astype(int)
    extra=draw-int(quota.sum());quota[np.argsort(-(desired-quota),kind='stable')[:extra]]+=1
    assert quota.sum()==draw and (quota<=sizes).all()
    return group,quota


def project(value,group,quota,draw):
    result=np.zeros(len(value));cap=1/draw
    for b,k in enumerate(quota):
        idx=np.flatnonzero(group==b)
        if k==0:continue
        if k==len(idx):result[idx]=cap;continue
        z=np.asarray(value)[idx];mass=k/draw;lo=z.min()-cap;hi=z.max()
        for _ in range(55):
            mid=(lo+hi)/2
            if np.clip(z-mid,0,cap).sum()>mass:lo=mid
            else:hi=mid
        result[idx]=np.clip(z-(lo+hi)/2,0,cap)
    return result


def uniform(group,quota,draw):
    result=np.zeros(len(group))
    for b,k in enumerate(quota):result[group==b]=k/draw/np.count_nonzero(group==b)
    return result


def pivotal(probabilities,rng):
    p=np.asarray(probabilities,float).copy();n=int(round(p.sum()))
    assert np.min(p)>=-1e-9 and np.max(p)<=1+1e-9 and abs(p.sum()-n)<1e-8
    p=np.clip(p,0,1);order=rng.permutation(len(p));active=None
    for j in order:
        if p[j]<1e-12:p[j]=0;continue
        if p[j]>1-1e-12:p[j]=1;continue
        if active is None:active=j;continue
        i=active;a=p[i];b=p[j];total=a+b
        if total<1:
            if rng.random()<a/total:p[i]=total;p[j]=0;active=i
            else:p[i]=0;p[j]=total;active=j
        else:
            if rng.random()<(1-b)/(2-total):p[i]=1;p[j]=total-1;active=j
            else:p[i]=total-1;p[j]=1;active=i
        if p[active]<1e-12:p[active]=0;active=None
        elif p[active]>1-1e-12:p[active]=1;active=None
    if active is not None:
        assert min(p[active],1-p[active])<1e-8
        p[active]=round(p[active])
    ids=np.flatnonzero(p>.5);assert len(ids)==n and len(np.unique(ids))==n
    return ids


def sample(weights,group,quota,draw,seed):
    rng=np.random.default_rng(seed);chosen=[]
    for b,k in enumerate(quota):
        idx=np.flatnonzero(group==b);local=pivotal(np.asarray(weights)[idx]*draw,rng)
        assert len(local)==k;chosen.extend(idx[local])
    ids=np.asarray(chosen,int);assert len(ids)==draw and len(np.unique(ids))==draw
    return ids[rng.permutation(draw)]


def edit_own(frame,own,seed):
    result=frame.copy()
    if own is not None:
        mask=result.indices==own
        result.data[mask]=np.random.default_rng(seed).binomial(result.data[mask],.1)
        result.eliminate_zeros()
    return result


def operators(x,target,axis):
    libraries=np.asarray(x.sum(1)).ravel().astype(float);scale=float(libraries.mean())
    derivative=np.zeros(len(target));derivative[axis]=50000/(1+50000*target[axis])
    def forward(w):return (np.asarray(x.T@w).ravel()-target*float(libraries@w))*derivative/scale
    def transpose(v):
        d=derivative*v
        return (np.asarray(x@d).ravel()-libraries*float(target@d))/scale
    return forward,transpose


def fit(frame,target,axis,draw,group=None,quota=None):
    libs=np.asarray(frame.sum(1)).ravel().astype(float)
    if group is None:group,quota=strata(libs,draw)
    u=uniform(group,quota,draw);p=u.copy();y=p.copy();n=len(p);naxis=int(np.count_nonzero(axis));assert naxis>0
    forward,transpose=operators(frame,target,axis)
    def value(w):
        z=forward(w)
        return .5*float(z@z)/naxis+.5*RIDGE*n*float((w-u)@(w-u))
    def gradient(w):return transpose(forward(w))/naxis+RIDGE*n*(w-u)
    initial=value(p);previous=initial;L=1.+RIDGE*n;t=1.;gap=float('inf');history=[]
    for iteration in range(250):
        gy=gradient(y);fy=value(y)
        for backtrack in range(40):
            trial=project(y-gy/L,group,quota,draw);d=trial-y;ftrial=value(trial)
            if ftrial<=fy+float(gy@d)+.5*L*float(d@d)+1e-15:break
            L*=2
        else:raise RuntimeError('Population balancing line search failed')
        if ftrial>previous+1e-14:
            y=p.copy();t=1.;continue
        old=p;p=trial;previous=ftrial;tn=(1+np.sqrt(1+4*t*t))/2;y=p+(t-1)/tn*(p-old);t=tn
        if iteration%10==0 or iteration==249:
            gp=gradient(p);vertex=np.zeros(n)
            for b,k in enumerate(quota):
                idx=np.flatnonzero(group==b);vertex[idx[np.argsort(gp[idx],kind='stable')[:k]]]=1/draw
            gap=float(gp@(p-vertex));history.append(dict(iteration=iteration+1,objective=previous,duality_gap=gap))
            if gap<=1e-8:break
    assert p.min()>=0 and p.max()<=1/draw+1e-12
    for b,k in enumerate(quota):assert abs(p[group==b].sum()-k/draw)<1e-10
    assert previous<=initial+1e-12
    weighted=np.asarray(frame.T@p).ravel();weighted/=weighted.sum()
    return p,dict(initial_objective=initial,final_objective=previous,duality_gap=gap,
                  converged=gap<=1e-8,iterations=iteration+1,effective_sample_size=float(1/(p@p)),
                  maximum_inclusion_probability=float(p.max()*draw),history=history,
                  weighted_bulk_profile_error=float(np.mean((np.log1p(50000*weighted[axis])-np.log1p(50000*target[axis]))**2)))


def permute_within(weights,group,seed=68068):
    result=weights.copy();rng=np.random.default_rng(seed)
    for b in np.unique(group):
        idx=np.flatnonzero(group==b);result[idx]=weights[idx[rng.permutation(len(idx))]]
    return result


def freeze():
    from run_v35_v36 import contexts
    A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
    deps=['v68_population_balancing.py','run_v68_experiment.py','test_v68_balancing.py',
          'v25_methods.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py','audit_transfer_validation.py','run_v35_v36.py']
    inputs=['experiments/v66_assay_strength/analysis/release_manifest.json',
            'experiments/v66_assay_strength/analysis/leaderboard_result.json',
            'experiments/v34_residual_flow/analysis/plan.json']
    for ctx in contexts():
        name=ctx['family'];inputs += ['experiments/v66_assay_strength/analysis/predictions/'+name+'.json',
                                     'experiments/v37_component_ensemble/analysis/cache/'+name+'/receipt.json']
    write(A/'plan.json',dict(version=68,dependencies={n:sha(ROOT/n) for n in deps},inputs={n:sha(ROOT/n) for n in inputs},
          primary='balanced_real',controls=['v63','uniform_real','shuffled_weights'],
          hypothesis='Replace synthetic transcript generation with optimized sampling of distinct real control cells. Preserve every non-own count exactly; only own target is independently thinned by binomial probability .1 before fitting/sampling. Predict responses as changes in population occupancy constrained to measured control states. No post-sampling moment refit or additional noise.',
          objective='Convex quadratic discrepancy in a linearized TS50k bulk profile: A_jg=(X_jg-library_j*target_probability_g)*[50000/(1+50000*target_probability_g)]/mean_library. Minimize ||sum_j p_j A_j||^2/(2*n_nonpanel_measured_genes) + .5*1e-4*N*||p-u||^2. Target probability is the frozen V63 desired bulk profile, including own factor once. Fit uses query controls only, never recipient treatment or reference outcomes.',
          constraints='Eight rank-depth strata defined on original query controls. Largest-remainder sample quotas proportional to original counts. Nonnegative weights, each <=1/draw, each stratum mass=quota/draw. Accelerated projected gradient with backtracking, max250iterations, exact convex dual-gap stopping1e-8. Sampling via randomized pivotal splitting within strata, no replacement, exact quotas. A capped weight implies ESS>=draw. Near-optimality is reported; no invalid or duplicate sample may be emitted.',
          frame='All fixed query-fit NTC cells in original left control partition, disjoint from real reference and all treatments. Existing five contexts and studies; source effects are the study-excluded frozen V63 predictions. No new independent study or blind-data claim. D3 unopened.',
          development='First two unique target cases per context, original generated-cell count (200 or400) matched to previous comparator, original seed; same treated populations and real reference across arms. V63 saved count outcomes may be reused only with matching hashes/truth. Ten target/context pairs; four arms. One400-cell zero-response null per context for unchanged V63 and uniform real-cell sampling. Counts are run even if convex-hull projection cannot match intended means, because the proposed gain concerns distributions.',
          guards='Null significant calls<=max(5,V63 null calls) in every context; paired equal-study profile error <=.99 V63; PDS loss<=.02, fidelity loss<=.005, reach loss<=.01, Jaccard loss<=.001, NMAE error no worse>1%; >=2 DE metrics improve; beats shuffled and uniform controls on profile and >=2 available DE metrics. Each study resolved; report undefined values rather than zero. No promotion or parameter tuning after outcome. Only passing development proceeds to remaining canonical target cases/two seeds and then independent reserve before packaging.',
          references=['https://doi.org/10.1080/01621459.2015.1023805','https://doi.org/10.1093/biomet/85.1.89'],
          scope='Independent balancing/sampling adaptation inspired by stable weighting and pivotal sampling, not a reproduction of their causal-inference claims. It cannot create expression states absent from the control convex hull.',
          D3_opened=False,automatic_upload=False))
    print('V68 protocol frozen',sha(A/'plan.json'),flush=True)


def check():
    p=read(A/'plan.json')
    for name,h in {**p['dependencies'],**p['inputs']}.items():assert sha(ROOT/name)==h,name
    return p


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze']);args=p.parse_args();globals()[args.stage]()
