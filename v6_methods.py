"""V6 count emission and unanchored source metrics, separated from response fitting."""
import numpy as np
from scipy import sparse, stats
import build_submission as build
from build_v5_submission import project_bulk

CANDIDATES = {
    'v4': {'blend': 0., 'calibration': 0., 'dual': False},
    'v4_dual': {'blend': 0., 'calibration': 1., 'dual': True},
    'v5_response_v4_generator': {'blend': 1., 'calibration': 0., 'dual': False},
    'v5_dual': {'blend': 1., 'calibration': 1., 'dual': True},
    'v4_cal10': {'blend': 0., 'calibration': .1, 'dual': False},
    'v4_cal25': {'blend': 0., 'calibration': .25, 'dual': False},
    'blend10_v4': {'blend': .1, 'calibration': 0., 'dual': False},
    'blend25_v4': {'blend': .25, 'calibration': 0., 'dual': False},
}


def moments(x):
    x = sparse.csr_matrix(x, dtype=np.float64)
    depth = np.asarray(x.sum(axis=1)).ravel()
    if (depth <= 0).any():
        raise ValueError('Empty cell')
    normalized = x.multiply((1/depth)[:,None]).tocsr()
    p = np.asarray(normalized.mean(axis=0)).ravel()
    b = np.asarray(x.sum(axis=0)).ravel()/depth.sum()
    variance = np.maximum(np.asarray(normalized.power(2).mean(axis=0)).ravel()-p*p,0)
    return depth, p, b, variance


def profiles(mean, bulk, effect, own, knockdown=.9):
    factor = np.exp2(np.clip(effect, -3, 3)); factor[own] = 1-knockdown
    p, b = mean*factor, bulk*factor
    return p/p.sum(), b/b.sum(), factor


def emit(x, rows, effect, own, seed, setting, atlas, requested_bulk=None, prepared=None):
    if prepared is None:
        depth, p0, b0, _ = moments(x)
        mean_raw = np.asarray(x.mean(axis=0)).ravel()
    else:
        depth, p0, b0, mean_raw = prepared
    p, b, factor = profiles(p0,b0,effect,own)
    if requested_bulk is not None:
        b = np.asarray(requested_bulk, dtype=float); b = b/b.sum()
    rng = np.random.default_rng(seed)
    initial = build.simulate_counts_zero_aware(sparse.csr_matrix(x), rows, factor, rng, mean_raw, depth.mean())
    if setting['calibration'] == 0:
        return initial, {'calibration': 0., 'requested_mean_L1': float(np.abs(moments(initial)[1]-p).sum()),
            'requested_bulk_L1': float(np.abs(moments(initial)[2]-b).sum()), 'projection_L1':0.}
    if setting['dual']:
        initial = sparse.csr_matrix(x)[rows].astype(np.int32)
    libs, pi, bi, _ = moments(initial)
    weight = setting['calibration']
    desired_p = (1-weight)*pi + weight*p
    desired_b = (1-weight)*bi + weight*b
    projected = project_bulk(desired_p,desired_b,libs)
    template = initial.toarray().astype(float)/libs[:,None]
    count = atlas.dual_moment_counts(template,desired_p,desired_b,depths=libs.astype(np.int64),seed=int(seed),iterations=200)
    np.testing.assert_array_equal(count.sum(axis=1),libs)
    _, pf, bf, _ = moments(count)
    diagnostic = {'calibration': weight,'requested_mean_L1':float(np.abs(pf-desired_p).sum()),
        'requested_bulk_L1':float(np.abs(bf-projected).sum()),'projection_L1':float(np.abs(projected-desired_b).sum())}
    if diagnostic['requested_mean_L1']>.03 or diagnostic['requested_bulk_L1']>.002:
        raise ValueError('V6 integer moment fit failed')
    return sparse.csr_matrix(count), diagnostic


def de(x, control, genes):
    """Actual two-sided asymptotic Wilcoxon, CPM >5 control filter, within-row BH."""
    x, control = sparse.csr_matrix(x), sparse.csr_matrix(control)
    xd = np.asarray(x.sum(axis=1)).ravel(); cd = np.asarray(control.sum(axis=1)).ravel()
    xn = x.multiply((1e6/xd)[:,None]).tocsr(); cn = control.multiply((1e6/cd)[:,None]).tocsr()
    xm, cm = np.asarray(xn.mean(axis=0)).ravel(), np.asarray(cn.mean(axis=0)).ravel()
    eligible = cm>5; columns = np.flatnonzero(eligible)
    p = np.ones(len(genes))
    for start in range(0,len(columns),256):
        col=columns[start:start+256]
        p[col]=stats.mannwhitneyu(xn[:,col].toarray(),cn[:,col].toarray(),axis=0,
            method='asymptotic',use_continuity=True).pvalue
    q = np.ones(len(genes)); q[eligible]=stats.false_discovery_control(p[eligible],method='bh')
    lfc = np.log2((xm+1e-9)/(cm+1e-9))
    return {'p':p,'q':q,'lfc':lfc,'eligible':eligible}


def direction_metrics(pred, truth, genes, own):
    valid = truth['eligible'] & (np.arange(len(genes))!=own)
    real = valid & (truth['q']<.05)
    adjudicable = valid & np.isfinite(truth['lfc']) & (truth['lfc']!=0)
    called = adjudicable & (pred['q']<.05)
    correct = called & (np.sign(pred['lfc'])==np.sign(truth['lfc']))
    budget = max(int(real.sum()),int(called.sum()))
    fidelity = float(correct.sum()/budget) if budget else np.nan
    union = real | (valid & (pred['q']<.05))
    jac = float((real & (pred['q']<.05)).sum()/union.sum()) if union.any() else 1.
    pool=np.flatnonzero(real & adjudicable)
    order=np.lexsort((np.asarray(genes)[pool],-np.abs(pred['lfc'][pool]),pred['p'][pool],pred['q'][pool],pred['q'][pool]>=.05))
    pool=pool[order]
    if len(pool):
        pure=np.cumsum(np.sign(pred['lfc'][pool])==np.sign(truth['lfc'][pool]))/np.arange(1,len(pool)+1)
        passes=np.flatnonzero(pure>=.9)
        reach=float((passes[-1]+1)/int(real.sum())) if len(passes) else 0.
    else:
        reach=np.nan
    nmae=float(np.abs(pred['lfc'][real]-truth['lfc'][real]).sum()/np.abs(truth['lfc'][real]).sum()) if real.any() else np.nan
    return {'fidelity':fidelity,'reach':reach,'jaccard':jac,'nmae':nmae,
        'predicted_calls':int(called.sum()),'reference_calls':int(real.sum())}


def profile_jackknife(x):
    x = sparse.csr_matrix(x).toarray().astype(float)
    depth=x.sum(axis=1); total=depth.sum(); counts=x.sum(axis=0)
    profile=np.log1p(50000*counts/total)
    if len(x)>1:
        deleted=np.log1p(50000*(counts[None,:]-x)/(total-depth)[:,None])
        jack=(len(x)-1)*deleted.var(axis=0)
    else:
        jack=np.zeros(x.shape[1])
    return profile,jack


def summarize_profiles(pred, truth, ctrl, pj, tj, cj, genes, targets, official_targets):
    """Unanchored PDS and panel corrected expression ratio on the available source panel."""
    own = np.array([int(np.flatnonzero(genes==t)[0]) for t in targets])
    mask = np.ones(pred.shape,dtype=bool); mask[np.arange(len(targets)),own]=False
    size = mask.sum(axis=1)
    means = np.divide((pred*mask).sum(axis=0),mask.sum(axis=0),out=np.zeros(pred.shape[1]),where=mask.sum(axis=0)>0)
    spread=((pred-means)**2*mask/size[:,None]).sum()
    predc=(pj*mask).sum(axis=1)/size; realc=(tj*mask).sum(axis=1)/size; ctrlc=(cj*mask).sum(axis=1)/size
    capped=np.minimum(predc,realc)
    rho=min(1.,float(spread/capped.sum())) if capped.sum()>0 else 1.
    numerator=(((pred-truth)**2*mask).sum(axis=1)/size-rho*capped-realc).sum()
    denominator=(((truth-ctrl)**2*mask).sum(axis=1)/size-realc-ctrlc).sum()
    expr=max(0.,float(numerator))/max(float(denominator),1e-12)
    columns=~np.isin(genes,official_targets)
    a,b=(pred-ctrl)[:,columns],(truth-ctrl)[:,columns]
    norm=np.linalg.norm(a,axis=1)[:,None]*np.linalg.norm(b,axis=1)[None,:]
    distance=1-np.divide(a@b.T,norm,out=np.zeros((len(targets),len(targets))),where=norm>0)
    diagonal=np.diag(distance)
    ranks=(distance<diagonal[:,None]).sum(axis=1)+.5*((distance==diagonal[:,None]).sum(axis=1)-1)
    pds=1-ranks/max(len(targets)-1,1)
    return {'expression_ratio':expr,'expression_numerator':float(numerator),'expression_denominator':float(denominator),
        'expression_rho':rho,'pds':float(pds.mean()) if len(targets)>1 else np.nan,'uncorrected_profile_mse':float(((pred-truth)**2*mask).sum()/mask.sum())}
