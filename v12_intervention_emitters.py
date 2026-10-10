"""Stage-two count mechanisms: response-linked variance and physical transcript updates."""
import numpy as np
from scipy import sparse
from v6_methods import moments
from v12_methods import emit as initial_emit

def emit(x,rows,effect,own,seed,emitter,prepared=None,atlas=None):
 if emitter in ['v11','gamma100','gamma25','empirical100']:return initial_emit(x,rows,effect,own,seed,emitter,prepared,atlas)
 if prepared is None:prepared=moments(x)
 depth,p,b,var=prepared[:4];rows=np.asarray(rows,int);libs=np.asarray(depth)[rows].astype(np.int64)
 factor=np.exp2(np.clip(effect,-3,3))
 if own is not None:factor[own]=.1
 wp=p*factor;wp/=wp.sum();wb=b*factor;wb/=wb.sum()
 rng=np.random.default_rng(seed);original=x[rows].toarray().astype(np.int64)
 if emitter=='thinning_immigration':
  # CRISPR response removes existing transcripts probabilistically, or adds
  # independent transcripts at the requested mean; zero-response genes keep
  # their sampled biological state and dropout. No expression pseudo-counts.
  initial=rng.binomial(original,np.minimum(factor,1)[None,:]).astype(float)
  increase=np.maximum(factor-1,0)*p
  initial+=rng.poisson(libs[:,None]*increase[None,:])
 elif emitter=='response_linked_variance':
  empirical=original*factor[None,:];empirical/=empirical.sum(axis=1,keepdims=True)
  expected=libs[:,None]*wp;compact=.35*expected+.65*rng.poisson(expected);compact/=compact.sum(axis=1,keepdims=True)
  # A single declared threshold: zero effects preserve empirical variability;
  # only appreciable predicted responses use the compact template.
  strength=np.clip(abs(effect)/.25,0,1)
  if own is not None:strength[own]=1
  initial=(1-strength[None,:])*empirical+strength[None,:]*compact
 else:raise ValueError(emitter)
 initial/=initial.sum(axis=1,keepdims=True)
 if atlas is None:
  from evaluate_atlas_transfer import load_atlas_model
  atlas,_=load_atlas_model()
 counts=atlas.dual_moment_counts(initial,wp,wb,depths=libs,seed=int(seed),iterations=200)
 np.testing.assert_array_equal(counts.sum(axis=1,dtype=np.int64),libs)
 assert np.issubdtype(counts.dtype,np.integer) and (counts>=0).all()
 rp=(counts/libs[:,None]).mean(axis=0);rb=counts.sum(axis=0,dtype=np.int64)/libs.sum()
 diag={'mean_probability_L1':float(abs(rp-wp).sum()),'bulk_probability_L1':float(abs(rb-wb).sum())}
 if diag['mean_probability_L1']>.03 or diag['bulk_probability_L1']>.031:raise ValueError('Moment constraints failed')
 return sparse.csr_matrix(counts),diag
