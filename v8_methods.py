"""V8 shared-control metrics, isolated response variants and calibrated integer emission."""
from pathlib import Path
import json
import numpy as np
from scipy import sparse
import build_submission as build
from audit_transfer_validation import pds, corrected_expression, profile_jackknife_chunked
from v6_methods import moments, de, direction_metrics
from v7_response_methods import collapse, combine, adaptive_shrink, BASE_WEIGHTS
from evaluate_atlas_transfer import load_atlas_model
from build_v5_submission import project_bulk
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v8_corrected_transfer';A=OUT/'analysis'
SEEDS=[20261008+8100,20261008+8119]
RESPONSE_GRID=[f'v4_a{a}' for a in [35,65,100]]+[f'{mode}_a{a}' for mode in ['uniform','uniform_eb'] for a in [25,50,75]]


def old_response(line=None,context=None):
    path=ROOT/'experiments/v4_lineage_response/analysis/source_components.npz'
    rec=json.loads(path.with_name('reference_model.json').read_text())
    if build.sha256(path)!=rec['source_components_sha256']:raise ValueError('V4 component hash mismatch')
    with np.load(path,allow_pickle=False) as z:
        if context=='A':return np.where(z['cd4_mask'][0],z['cd4'][0],np.where(z['k562_mask'],z['k562'],0))
        condition=1 if context=='B' else 2
        transferred=.5*z['cd4_specific'][condition]+z['k562']
        if line=='H1':return transferred
        return np.where(z['h1_mask'],z['h1'],transferred)+z['h1_prior']


def response_bank(data,excluded=None,context=None):
    train=np.ones(len(data['truth']),bool) if excluded is None else data['row_families']!=excluded
    cube=collapse(data,train)
    mean,var,valid=combine(cube,BASE_WEIGHTS)
    eb,mixture=adaptive_shrink(mean,var,valid)
    old=old_response(line=excluded,context=context)
    bank={}
    for name in RESPONSE_GRID:
        mode,amplitude=name.rsplit('_a',1)
        value=old if mode=='v4' else mean if mode=='uniform' else eb
        bank[name]=value*(int(amplitude)/100)
    return bank,{'excluded_family':excluded,'train_rows':np.flatnonzero(train).tolist(),
                 'train_families':sorted(set(data['row_families'][train].tolist())),
                 'pooling_weights':BASE_WEIGHTS.tolist(),'uniform_eb_prior':mixture.tolist()}


def gated(effect,query,own):
    gate=np.clip(query/5,0,1)
    result=np.clip(effect,-3,3)*gate*gate[own]
    result[own]=0
    return result


def emit(x,rows,effect,own,seed,generator='legacy',prepared=None,atlas=None):
    """Same response; vary only count realization. None means a true control null."""
    if prepared is None:
        depth,p,b,_=moments(x);raw=np.asarray(x.mean(axis=0)).ravel()
    else:depth,p,b,raw=prepared
    factor=np.exp2(np.clip(effect,-3,3))
    if own is not None:factor[own]=.1
    if generator=='legacy':
        counts=build.simulate_counts_zero_aware(x,rows,factor,np.random.default_rng(seed),raw,depth.mean())
        return counts,{'generator':generator}
    if generator not in ['dual','dual_half']:raise ValueError('Unknown V8 generator')
    template=sparse.csr_matrix(x)[rows]
    libs=np.asarray(template.sum(axis=1)).ravel().astype(np.int64)
    wanted_p=p*factor;wanted_p/=wanted_p.sum()
    wanted_b=b*factor;wanted_b/=wanted_b.sum()
    projected=project_bulk(wanted_p,wanted_b,libs)
    projection=float(np.abs(projected-wanted_b).sum())
    if projection>.03:raise ValueError('Infeasible mean/bulk constraints')
    probabilities=template.toarray().astype(float)/libs[:,None]
    if generator=='dual_half':probabilities=.5*probabilities+.5*wanted_p
    if atlas is None:atlas,_=load_atlas_model()
    counts=atlas.dual_moment_counts(probabilities,wanted_p,wanted_b,depths=libs,seed=int(seed),iterations=200)
    np.testing.assert_array_equal(counts.sum(axis=1,dtype=np.int64),libs)
    realized_p=np.mean(counts/libs[:,None],axis=0)
    realized_b=counts.sum(axis=0,dtype=np.int64)/libs.sum()
    diag={'generator':generator,'mean_probability_L1':float(np.abs(realized_p-wanted_p).sum()),
          'bulk_probability_L1':float(np.abs(realized_b-wanted_b).sum()),
          'projected_bulk_probability_L1':float(np.abs(realized_b-projected).sum()),
          'bulk_projection_L1':projection}
    if diag['mean_probability_L1']>.03 or diag['projected_bulk_probability_L1']>.002:
        raise ValueError('Integer moment calibration failed')
    return sparse.csr_matrix(counts),diag


def panel_metrics(pred,truth,control,pj,tj,cj,genes,target_names,panel_targets,reference_control=None):
    own=np.array([np.flatnonzero(genes==t)[0] for t in target_names])
    mask=~np.isin(genes,panel_targets)
    if reference_control is None:
        discrimination=pds(pred[:,mask],truth[:,mask],control[mask])
    else:
        # Independent controls only for replicate diagnostics, not model scoring.
        from audit_transfer_validation import distance_rank
        discrimination=distance_rank((pred-control)[:,mask],(truth-reference_control)[:,mask])
    result=corrected_expression(pred,truth,control if reference_control is None else reference_control,pj,tj,cj,own)
    return {'pds':discrimination,**result}
