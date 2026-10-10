"""V12 isolated count distributions and corrected competition-style CPU metrics."""
import numpy as np
from scipy import sparse,stats
from v6_methods import moments
from v6_methods import direction_metrics as legacy_direction
from v11_methods import emit as compact_emit,SETTINGS

def de(x,control,genes):
 x=sparse.csr_matrix(x);control=sparse.csr_matrix(control)
 xd=np.asarray(x.sum(axis=1)).ravel();cd=np.asarray(control.sum(axis=1)).ravel()
 assert (xd>0).all() and (cd>0).all()
 xn=x.multiply((1e6/xd)[:,None]).tocsr();cn=control.multiply((1e6/cd)[:,None]).tocsr()
 xm=np.asarray(xn.mean(axis=0)).ravel();cm=np.asarray(cn.mean(axis=0)).ravel()
 eligible=cm>5;cols=np.flatnonzero(eligible);p=np.ones(len(genes));q=p.copy()
 for start in range(0,len(cols),256):
  c=cols[start:start+256]
  p[c]=stats.mannwhitneyu(xn[:,c].toarray(),cn[:,c].toarray(),axis=0,method='asymptotic',use_continuity=False).pvalue
 p=np.nan_to_num(p,nan=1.);q[eligible]=stats.false_discovery_control(p[eligible],method='bh')
 return {'p':p,'q':q,'lfc':np.log2((xm+1e-9)/(cm+1e-9)),'eligible':eligible}

def direction_metrics(pred,truth,genes,own):
 result=legacy_direction(pred,truth,genes,own)
 real=truth['eligible']&(np.arange(len(genes))!=own)&(truth['q']<.05)&np.isfinite(truth['lfc'])
 denom=abs(truth['lfc'][real]).sum()
 result['nmae']=float(abs(np.nan_to_num(pred['lfc'][real],nan=0.,posinf=0.,neginf=0.)-truth['lfc'][real]).sum()/denom) if real.sum()>=10 and np.isfinite(denom) and denom>0 else np.nan
 if result['reference_calls']>0 and not np.isfinite(result['reach']):result['reach']=0.
 result['nmae_eligible']=bool(np.isfinite(result['nmae']))
 return result

def emit(x,rows,effect,own,seed,emitter,prepared=None,atlas=None):
 if prepared is None:prepared=moments(x)
 depth,p,b,var=prepared[:4];rows=np.asarray(rows,dtype=int);libs=np.asarray(depth)[rows].astype(np.int64)
 if emitter=='v11':return compact_emit(x,rows,effect,own,seed,SETTINGS['consensus20_noise65'],prepared=(depth,p,b),atlas=atlas)
 factor=np.exp2(np.clip(effect,-3,3))
 if own is not None:factor[own]=.1
 wp=p*factor;wp/=wp.sum();wb=b*factor;wb/=wb.sum()
 rng=np.random.default_rng(seed)
 if emitter.startswith('gamma'):
  # Fit biological composition dispersion from the complete FIT population,
  # never the sampled template depths or the held-out reference.
  excess=np.maximum(var-p*(prepared[4] if len(prepared)>4 else np.mean(1/np.asarray(depth))),0)
  alpha=np.clip(np.divide(p*p,excess,out=np.full(len(p),1e8),where=excess>0),1e-3,1e8)
  latent=rng.gamma(alpha[None,:],(wp/alpha)[None,:],size=(len(rows),len(p)))
  initial=rng.poisson(libs[:,None]*latent).astype(float)
  if emitter=='gamma25':
   initial/=np.maximum(initial.sum(axis=1,keepdims=True),1)
   expected=libs[:,None]*wp
   compact=.35*expected+.65*rng.poisson(expected);compact/=compact.sum(axis=1,keepdims=True)
   initial=.25*initial+.75*compact
 elif emitter=='empirical100':
  initial=x[rows].toarray().astype(float)*factor[None,:]
 else:raise ValueError(emitter)
 initial/=np.maximum(initial.sum(axis=1,keepdims=True),1)
 if atlas is None:
  from evaluate_atlas_transfer import load_atlas_model
  atlas,_=load_atlas_model()
 counts=atlas.dual_moment_counts(initial,wp,wb,depths=libs,seed=int(seed),iterations=200)
 np.testing.assert_array_equal(counts.sum(axis=1,dtype=np.int64),libs)
 assert np.issubdtype(counts.dtype,np.integer) and (counts>=0).all()
 rp=np.mean(counts/libs[:,None],axis=0);rb=counts.sum(axis=0,dtype=np.int64)/libs.sum()
 diag={'mean_probability_L1':float(abs(rp-wp).sum()),'bulk_probability_L1':float(abs(rb-wb).sum())}
 if diag['mean_probability_L1']>.03 or diag['bulk_probability_L1']>.031:raise ValueError('Moment constraints failed')
 return sparse.csr_matrix(counts),diag
