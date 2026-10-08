"""V9 profile calibration initialized from zero-aware simulated perturbations."""
from pathlib import Path
import numpy as np
from scipy import sparse
from v6_methods import moments
from evaluate_atlas_transfer import load_atlas_model
from build_v5_submission import project_bulk
import build_submission as build
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'experiments/v9_activation_preserving';A=OUT/'analysis'
SEEDS=[20269108,20269127]
GENERATORS=['legacy','legacy_fit25','legacy_fit50','legacy_fit100']


def emit(x,rows,effect,own,seed,generator='legacy',prepared=None,atlas=None):
    if prepared is None:
        depth,p,b,_=moments(x);raw=np.asarray(x.mean(axis=0)).ravel()
    else:depth,p,b,raw=prepared
    factor=np.exp2(np.clip(effect,-3,3))
    if own is not None:factor[own]=.1
    initial=build.simulate_counts_zero_aware(x,rows,factor,np.random.default_rng(seed),raw,depth.mean())
    if generator=='legacy':return initial,{'generator':generator}
    if generator not in GENERATORS:raise ValueError('Unknown V9 generator')
    strength=int(generator.replace('legacy_fit',''))/100
    libs,pi,bi,_=moments(initial);libs=libs.astype(np.int64)
    wanted_p=p*factor;wanted_p/=wanted_p.sum();wanted_b=b*factor;wanted_b/=wanted_b.sum()
    desired_p=(1-strength)*pi+strength*wanted_p
    desired_b=(1-strength)*bi+strength*wanted_b
    projected=project_bulk(desired_p,desired_b,libs)
    if np.abs(projected-desired_b).sum()>.03:raise ValueError('Infeasible V9 moment request')
    probabilities=initial.toarray().astype(float)/libs[:,None]
    if atlas is None:atlas,_=load_atlas_model()
    counts=atlas.dual_moment_counts(probabilities,desired_p,desired_b,depths=libs,seed=int(seed),iterations=200)
    np.testing.assert_array_equal(counts.sum(axis=1,dtype=np.int64),libs)
    pf=counts/libs[:,None];realized_p=pf.mean(axis=0);realized_b=counts.sum(axis=0,dtype=np.int64)/libs.sum()
    diag={'generator':generator,'strength':strength,'mean_probability_L1':float(abs(realized_p-desired_p).sum()),
          'bulk_probability_L1':float(abs(realized_b-desired_b).sum()),
          'projected_bulk_probability_L1':float(abs(realized_b-projected).sum()),
          'bulk_projection_L1':float(abs(projected-desired_b).sum()),
          'mean_detected_genes_before_calibration':float(initial.getnnz(axis=1).mean()),
          'mean_detected_genes_after_calibration':float(np.count_nonzero(counts,axis=1).mean())}
    if diag['mean_probability_L1']>.03 or diag['projected_bulk_probability_L1']>.002:raise ValueError('V9 integer moment fit failed')
    return sparse.csr_matrix(counts),diag
