"""V69 independent Gamma-Poisson count bridge and source dispersion transfer.

No code copied from external implementations. Zero response is count identity.
"""
from pathlib import Path
import numpy as np
from scipy import sparse
from v57_pie import read,write,sha
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v69_count_bridge';A=OUT/'analysis'
SOURCES=('H1','HCT116','HEK293T','CD4_Rest_D2','CD4_Stim48_D2')
ARMS=('mean_bridge','shape_bridge','hybrid25')
ALPHA_MIN=.005;ALPHA_MAX=20.


def study(name):
    return 'CD4' if name.startswith('CD4') else ('Orion' if name in ('HCT116','HEK293T') else name)


def parameters(x):
    """Depth-conditioned method of moments; a working NB model, not an MLE."""
    x=sparse.csr_matrix(x,dtype=np.float64);n=x.shape[0]
    depth=np.asarray(x.sum(1)).ravel();assert n>=2 and (depth>0).all()
    total=np.asarray(x.sum(0)).ravel();p=total/depth.sum();d2=float(depth@depth)
    xx=np.asarray(x.multiply(x).sum(0)).ravel();dx=np.asarray(x.T@depth).ravel()
    denom=p*p*d2
    numerator=xx-total-2*p*dx+denom
    alpha=np.divide(numerator,denom,out=np.zeros_like(p),where=denom>0)
    return dict(p=p,alpha=np.clip(alpha,ALPHA_MIN,ALPHA_MAX),n=n,
                total=total,mean_library=float(depth.mean()),depth=depth)


def bridge(template,effect,own,seed,params,delta=None,strength=.5):
    """Poisson coupling via posterior Gamma intensities; allows varying libraries.

    Under the working model, X|lambda~Poisson(lambda), lambda~Gamma(r0,mu0/r0).
    Draw lambda|X; beta/gamma shape coupling produces Gamma(r1,mu1/r1).
    Thin X when lambda decreases; superpose independent counts when it increases.
    """
    original=sparse.csr_matrix(template).toarray().astype(np.int64)
    assert (original>=0).all()
    effect=np.asarray(effect,float);alpha=np.asarray(params['alpha'],float)
    change=np.zeros_like(alpha) if delta is None else strength*np.asarray(delta,float)
    if own is None and np.count_nonzero(effect)==0 and np.count_nonzero(change)==0:
        return sparse.csr_matrix(original),dict(zero_response_identity=True,changed_fraction=0.)
    factor=np.exp2(np.clip(effect,-3,3))
    if own is not None:factor[own]=.1;change=change.copy();change[own]=0.
    p=np.asarray(params['p'],float);closure=float(p@factor);assert closure>0
    factor=factor/closure
    r0=1/alpha;r1=1/np.clip(alpha*np.exp(np.clip(change,-np.log(4),np.log(4))),ALPHA_MIN,ALPHA_MAX)
    rng=np.random.default_rng(seed);result=original.copy();libraries=original.sum(1)
    for start in range(0,len(p),512):
        end=min(start+512,len(p));valid=p[start:end]>0
        if not valid.any():continue
        ix=np.arange(start,end)[valid];obs=original[:,ix];mu=libraries[:,None]*p[ix]
        a=r0[ix][None,:];b=r1[ix][None,:];ratio=factor[ix][None,:]*a/b
        latent=rng.gamma(a+obs,mu/(a+mu));shifted=latent.copy()
        higher=(r1[ix]-r0[ix])>1e-10;lower=(r0[ix]-r1[ix])>1e-10
        if higher.any():
            shifted[:,higher]+=rng.gamma((b-a)[:,higher],mu[:,higher]/a[:,higher])
        if lower.any():
            shifted[:,lower]*=rng.beta(b[:,lower],(a-b)[:,lower],size=(len(obs),int(lower.sum())))
        shifted*=ratio
        keep=np.minimum(np.divide(shifted,latent,out=np.ones_like(latent),where=latent>0),1.)
        removed=rng.binomial(obs,np.clip(1-keep,0,1));added=rng.poisson(np.maximum(shifted-latent,0))
        result[:,ix]=obs-removed+added
    assert (result>=0).all() and (result.sum(1)>0).all()
    return sparse.csr_matrix(result),dict(zero_response_identity=False,
        changed_fraction=float(abs(result-original).sum()/max(original.sum(),1)),
        expected_library_closure=closure,mean_library_ratio=float(result.sum()/original.sum()))


def check():
    plan=read(A/'plan.json')
    for name,h in plan['dependencies'].items():assert sha(ROOT/name)==h,name
    for name,h in plan['inputs'].items():assert sha(ROOT/name)==h,name
    return plan


def fit_source(name,genes,targets):
    from v25_methods import raw
    path=A/'source_shapes'/(name+'.npz');path.parent.mkdir(exist_ok=True)
    if path.exists():
        rec=read(path.with_suffix('.json'));assert sha(path)==rec['sha256'];return
    x,meta,measured,source=raw(name,genes);label=meta.gene_target.to_numpy(str)
    controls=np.flatnonzero(label=='Non-Targeting');panel=np.isin(genes,targets)
    shifts=np.zeros((len(targets),len(genes)),np.float32);weight=np.zeros_like(shifts);rows=[]
    batch=meta['sample'].astype(str).to_numpy() if 'sample' in meta else np.full(len(meta),'unknown')
    for ti,target in enumerate(targets):
        ids=np.flatnonzero(label==target)
        allowed=[b for b in np.unique(batch[ids]) if np.count_nonzero(batch[controls]==b)>=20]
        ids=ids[np.isin(batch[ids],allowed)];ctrl=controls[np.isin(batch[controls],allowed)]
        if len(ids)<80 or len(ctrl)<80:continue
        # Deterministic alternating cells within sorted source rows; sampling views,
        # not guide-disjoint or independent biological replicates.
        cp=parameters(x[ctrl]);fits=[parameters(x[ids[k::2]]) for k in (0,1)]
        dd=[np.clip(np.log(f['alpha']/cp['alpha']),-np.log(4),np.log(4)) for f in fits]
        agree=np.sign(dd[0])==np.sign(dd[1]);delta=(dd[0]+dd[1])/2
        eligible=measured&~panel&(cp['total']>=80)&(fits[0]['total']>=40)&(fits[1]['total']>=40)
        confidence=(len(ids)/(len(ids)+200))*np.minimum(1,cp['total']/400)*agree*eligible
        shifts[ti]=delta;weight[ti]=confidence
        rows.append(dict(target=str(target),treated_cells=len(ids),matched_controls=len(ctrl),
                         supported_coordinates=int(np.count_nonzero(confidence))))
    np.savez_compressed(path,delta=shifts,weight=weight,genes=np.asarray(genes,str),targets=np.asarray(targets,str))
    write(path.with_suffix('.json'),dict(sha256=sha(path),source=name,study=study(name),source_record=source,
          targets=rows,plan_sha256=sha(A/'plan.json'),split='Alternating source rows; shared-guide sampling views'))
    print('V69 fitted dispersion source',name,len(rows),flush=True)


def pooled(excluded=(),shuffle=False):
    from v62_conserved_response import data
    d=data();g=d['genes'];t=d['targets'];num={};den={};used=[]
    for name in SOURCES:
        group=study(name)
        if group in excluded:continue
        path=A/'source_shapes'/(name+'.npz');rec=read(path.with_suffix('.json'));assert sha(path)==rec['sha256']
        with np.load(path,allow_pickle=False) as z:
            np.testing.assert_array_equal(z['genes'],g);np.testing.assert_array_equal(z['targets'],t)
            delta=z['delta'].astype(float);weight=z['weight'].astype(float)
        num[group]=num.get(group,0)+delta*weight;den[group]=den.get(group,0)+weight;used.append(name)
    numerator=np.zeros((len(t),len(g)));denominator=np.zeros_like(numerator)
    for group in num:
        # Each study has at most one vote; a study with two lines is not doubled.
        w=np.minimum(den[group],1.);value=np.divide(num[group],den[group],out=np.zeros_like(num[group]),where=den[group]>0)
        numerator+=w*value;denominator+=w
    result=np.divide(numerator,denominator+1.,out=np.zeros_like(numerator),where=denominator>0)
    if shuffle:result=np.roll(result,1,axis=0)
    return result.astype(np.float32),dict(excluded_studies=list(excluded),training_contexts=used,
                                         shrinkage='One zero-shift prior study; study-balanced reliability weights')
