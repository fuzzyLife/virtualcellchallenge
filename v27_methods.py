"""V27 control-anchored neural response changes; historical arrays stay frozen."""
from pathlib import Path
import json
import numpy as np
from scipy import sparse
from v26_methods import (ROOT,geometry,decode,descriptor,fit_probabilities,integerize,
                         preserve_depth_margins,repair_cpm,load_decoder as old_decoder)
from v22_methods import write,digest,logcpm
import build_submission as b
OUT=ROOT/'experiments/v27_control_anchored'; A=OUT/'analysis'
ARMS=['v24','v26','fallback_only','anchor_only','calibrated_anchor','relative','calibrated_relative','panel_relative','target_shuffled']
FINALIST='panel_relative'

def load_plan():
    p=json.loads((A/'plan.json').read_text())
    for n,h in p['dependencies'].items():assert b.sha256(ROOT/n)==h,'Frozen V27 dependency changed: '+n
    for n,h in p['inputs'].items():assert b.sha256(ROOT/n)==h,'Frozen V27 input changed: '+n
    return p

def load_model(key,panel=False):
    kind='panel' if panel else 'conditional'
    f=A/'models'/f'{key}_{kind}.npz'
    if not f.exists() and not panel:return old_decoder(key,'conditional')
    rec=json.loads(f.with_suffix('.json').read_text())
    assert b.sha256(f)==rec['sha256'] and rec['plan_sha256']==b.sha256(A/'plan.json')
    with np.load(f,allow_pickle=False) as z:
        return {'weights':{n[3:]:z[n].copy() for n in z.files if n.startswith('w__')},
                **{n:z[n].copy() for n in ['output_mask','state_axis','state_basis','state_scale']},
                'geometry':geometry(key),'record':rec}

def relative_delta(state,effect,p,model,seed,override=None):
    # Common latent draws isolate the change in target condition. They do not
    # establish counterfactual pairing between source and query cells.
    z=np.random.default_rng(np.random.SeedSequence([seed,27027])).normal(size=(len(state),16)).astype(np.float32)
    c=descriptor(effect,p,model['geometry']) if override is None else override
    c=np.repeat(c[None,:],len(state),axis=0)
    delta=np.clip(decode(model,z,state,c)-decode(model,z,state,np.zeros_like(c)),-1,1)
    delta[:,~model['output_mask']]=0
    return delta

def control_calibration(template,fit,p,seed):
    # Only disjoint fitting NTCs enter this empirical quantile map. Random tie
    # breaking avoids identical row-order assignments for all zero-count genes.
    x=template.toarray().astype(float); libs=x.sum(1); fl=np.asarray(fit.sum(1)).ravel()
    v=fit.toarray().astype(float)/fl[:,None]; q=x/libs[:,None]
    rng=np.random.default_rng(seed+27027); order=np.lexsort((rng.uniform(size=q.shape),q),axis=0)
    quant=(np.arange(len(q))+.5)/len(q)
    values=np.quantile(v,quant,axis=0,method='inverted_cdf')
    result=np.empty_like(q); np.put_along_axis(result,order,values,axis=0)
    mean=result.mean(0); result*=np.divide(p,mean,out=np.ones_like(p),where=mean>0)
    result/=np.maximum(result.sum(1,keepdims=True),1e-15)
    return result

def anchored_initial(template,effect,own,p,delta=None,calibrated=None):
    libs=np.asarray(template.sum(1)).ravel()
    q=template.toarray().astype(float)/libs[:,None] if calibrated is None else calibrated.copy()
    factor=np.exp2(np.clip(effect,-3,3))
    if own is not None:factor[own]=.1
    q*=factor[None,:]
    if delta is not None:q*=np.exp(delta)
    return q/np.maximum(q.sum(1,keepdims=True),1e-15)

def project(initial,reference,seed):
    # Same margins as V26, without re-sampling observed measurement counts.
    ref=reference.toarray().astype(np.int32) if sparse.issparse(reference) else np.asarray(reference,np.int32)
    libs=ref.sum(1,dtype=np.int64); rp=(ref/libs[:,None]).mean(0); rb=ref.sum(0)/libs.sum()
    prob=fit_probabilities(initial,libs,rp,rb)
    x,n=integerize(prob,libs,ref.sum(0,dtype=np.int64),seed)
    x,m=preserve_depth_margins(x,ref,seed+1)
    x,k=repair_cpm(x,ref)
    drift=float(abs((x/libs[:,None]).mean(0)-rp).sum())
    return sparse.csr_matrix(x),{'mean_reference_L1':drift,'bulk_exact':True,
             'rounding_repair_molecules':n,'depth_repair_molecules':m,'CPM_repair_molecules':k}

def fallback_only(template,effect,own,p,seed,reference,model,state):
    # V26 learned logits and random measurement layer; fallback stays exactly at
    # its V24 counts through per-depth-stratum restricted integer repairs.
    from v26_methods import learned_initial
    ref=reference.toarray().astype(np.int32); libs=ref.sum(1); rp=(ref/libs[:,None]).mean(0)
    delta=learned_initial(template,state,effect,p,model,seed)
    q=rp[None,:]*np.exp(delta); active=model['output_mask'].copy()
    if own is not None:active[own]=False
    q[:,~active]=ref[:,~active]/libs[:,None]
    q/=q.sum(1,keepdims=True)
    q=fit_probabilities(q,libs,rp,ref.sum(0)/libs.sum())
    fixed=ref.copy(); remaining=libs-ref[:,~active].sum(1)
    qp=q[:,active]; qp/=np.maximum(qp.sum(1,keepdims=True),1e-15)
    rng=np.random.default_rng(seed+260260)
    sample=np.array([rng.multinomial(int(n),v) for n,v in zip(remaining,qp)],dtype=np.int32)
    moved=0
    for j,rows in enumerate(np.array_split(np.argsort(libs,kind='stable'),16)):
        qpart=sample[rows]/np.maximum(remaining[rows,None],1)
        part,n=integerize(qpart,remaining[rows],ref[rows][:,active].sum(0),seed+j)
        fixed[np.ix_(rows,np.flatnonzero(active))]=part; moved+=n
    assert np.array_equal(fixed[:,~active],ref[:,~active])
    np.testing.assert_array_equal(fixed.sum(1),libs); np.testing.assert_array_equal(fixed.sum(0),ref.sum(0))
    drift=float(abs((fixed/libs[:,None]).mean(0)-rp).sum())
    return sparse.csr_matrix(fixed),{'mean_reference_L1':drift,'bulk_exact':True,'fallback_exact':True,'depth_repair_molecules':moved}

def emit(template,effect,own,p,seed,reference,arm,model,state,fit=None,override=None):
    if arm=='fallback_only':return fallback_only(template,effect,own,p,seed,reference,model,state)
    calibrated=control_calibration(template,fit,p,seed) if arm in ['calibrated_anchor','calibrated_relative','panel_relative','target_shuffled'] else None
    delta=None if arm in ['anchor_only','calibrated_anchor'] else relative_delta(state,effect,p,model,seed,override)
    if own is not None and delta is not None:delta[:,own]=0
    q=anchored_initial(template,effect,own,p,delta,calibrated)
    return project(q,reference,seed+270270)
