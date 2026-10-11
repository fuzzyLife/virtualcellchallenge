"""V16 cross-family response shrinkage and zero-response identity count transport."""
from pathlib import Path
import json,gc
import numpy as np
from scipy import sparse,stats
from v15_ensemble_methods import cube,pool
from v7_response_methods import FAMILIES,bh
from v8_methods import gated
from v6_methods import moments
from v11_methods import emit as compact,SETTINGS
import build_submission as build
ROOT=build.ROOT;OUT=ROOT/'experiments/v16_control_preserving';A=OUT/'analysis'
CANDIDATES={'v13':('v13','compact'),'calibrated_compact':('calibrated','compact'),'balanced':('v13','balanced'),'calibrated_balanced':('calibrated','balanced'),'residual10':('v13','residual10'),'calibrated_residual10':('calibrated','residual10'),'uniform_compact':('matched','compact'),'uniform_balanced':('matched','balanced'),'uniform_residual10':('matched','residual10')}

def source_arrays(data,excluded):
    y,var,valid,own=cube(data,excluded);base=pool(y,var,valid,np.array([2.,.25,2.,1.,1.]),own)
    positive=np.zeros_like(base,dtype=np.uint8);same=np.zeros_like(base,dtype=float);total=valid.sum(0)
    for fi in range(len(FAMILIES)):
        for ti in range(len(own)):
            mask=valid[fi,ti].copy();mask[own[ti]]=False
            p=2*stats.norm.sf(abs(y[fi,ti])/np.sqrt(np.maximum(var[fi,ti],.0025)));q=bh(p,mask)
            agree=np.sign(y[fi,ti])==np.sign(base[ti]);positive[ti]+=(mask&(q<.05)&agree)
            same[ti]+=mask&agree
    agreement=np.divide(same,total,out=np.zeros_like(same),where=total>0)
    return {'base':base,'votes':positive,'agreement':agreement,'own':own}

def source_cache(data,excluded):
    folder=A/'source_banks_cache';folder.mkdir(parents=True,exist_ok=True);key='_'.join(sorted(excluded)) or 'all';file=folder/(key+'.npz');side=file.with_suffix('.json')
    plan=json.loads((A/'plan.json').read_text());inputs={'source':plan['source_cache_sha256'],'code':build.sha256(Path(__file__)),'excluded':sorted(excluded)}
    if side.exists():
        rec=json.loads(side.read_text());assert rec['inputs']==inputs and rec['sha256']==build.sha256(file)
        with np.load(file,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
    result=source_arrays(data,excluded)
    if len(excluded)==1:
        prior=ROOT/'experiments/v14_context_interactions/analysis/models'/(excluded[0]+'.npz')
        with np.load(prior,allow_pickle=False) as z:np.testing.assert_array_equal(result['base'],z['base'])
    np.savez_compressed(file,**result);build.write_json(side,{'inputs':inputs,'sha256':build.sha256(file),'training_families':[str(f) for f in FAMILIES if f not in excluded]})
    print('V16 source bank fitted:',key,flush=True);return result

def bins(effect,bank):
    magnitude=np.digitize(abs(effect),[.05,.2]);return (np.minimum(bank['votes'],2)*6+magnitude*2+(bank['agreement']>.75)).astype(np.uint8)

def fit_calibration(data,excluded,reserved=(),provider=source_cache):
    numerator=np.zeros(18);denominator=np.zeros(18);counts=np.zeros(18,dtype=np.int64);folds=[]
    nonpanel=~np.isin(data['genes'],data['targets'])
    for family in FAMILIES:
        if family in excluded:continue
        inner=sorted(set(excluded)|{str(family)});bank=provider(data,inner)
        rows=np.flatnonzero((data['row_families']==family)&~np.isin(data['targets'][data['row_targets']],reserved));names=np.unique(data['row_targets'][rows]);n=np.zeros(18);d=np.zeros(18);c=np.zeros(18,dtype=np.int64)
        for ti in names:
            selected=rows[data['row_targets'][rows]==ti]
            for ri in selected:
                query=np.nan_to_num(data['profiles'][data['row_contexts'][ri]],nan=0.);pred=gated(bank['base'][ti].copy(),query,int(bank['own'][ti]));idx=bins(pred,{'votes':bank['votes'][ti],'agreement':bank['agreement'][ti]})
                valid=data['valid'][ri]&nonpanel&np.isfinite(data['truth'][ri]);weight=1/max(int(valid.sum()),1)/len(selected)/max(len(names),1)
                n+=np.bincount(idx[valid],weights=(pred[valid]*data['truth'][ri,valid])*weight,minlength=18)
                d+=np.bincount(idx[valid],weights=pred[valid]**2*weight,minlength=18);c+=np.bincount(idx[valid],minlength=18)
        numerator+=n;denominator+=d;counts+=c;folds.append({'family':str(family),'excluded_predictor_families':inner,'target_count':len(names),'numerator':n.tolist(),'denominator':d.tolist()});del bank;gc.collect()
    ridge=.25*float(denominator.sum())/18;gain=np.divide(numerator+ridge,denominator+ridge,out=np.ones(18),where=denominator+ridge>0);gain=np.clip(gain,0,1)
    # At least two independently significant source families and >=0.2 log2
    # retain full V13 response; the calibration only suppresses weaker evidence.
    gain[16:18]=1.
    return {'excluded':list(excluded),'training_families':[str(f) for f in FAMILIES if f not in excluded],'reserved_label_targets_excluded':list(reserved),'gains':gain.tolist(),'ridge':ridge,'numerator':numerator.tolist(),'denominator':denominator.tolist(),'observations_per_bin':counts.tolist(),'inner_folds':folds,'V13_anchor_weight':.5}

def predict(bank,model,query):
    raw=np.asarray([gated(v.copy(),np.nan_to_num(query,nan=0.),int(o)) for v,o in zip(bank['base'],bank['own'])],dtype=np.float32)
    idx=bins(raw,bank);factor=.5+.5*np.asarray(model['gains'])[idx];cal=(raw*factor).astype(np.float32)
    scale=np.divide(abs(cal).sum(1),abs(raw).sum(1),out=np.ones(len(raw)),where=abs(raw).sum(1)>0)
    matched=(raw*scale[:,None]).astype(np.float32)
    return {'v13':raw,'calibrated':cal,'matched':matched},factor

def balanced_counts(template,factor,seed):
    original=sparse.csr_matrix(template).toarray().astype(np.int64);libs=original.sum(1);expected=original*factor[None,:];expected*=np.divide(libs,expected.sum(1))[:,None]
    rate=np.divide(expected,original,out=np.ones_like(expected),where=original>0)
    rng=np.random.default_rng(seed);removed=rng.binomial(original,np.maximum(1-np.minimum(rate,1),0));count=original-removed
    extra=np.maximum(expected-original,0)
    for i,total in enumerate(removed.sum(1)):
        if total:
            probability=extra[i]/extra[i].sum();count[i]+=rng.multinomial(int(total),probability)
    np.testing.assert_array_equal(count.sum(1),libs);assert (count>=0).all()
    return count,expected

def residual_initial(expected,poisson,empirical,alpha=.1):
    if not 0<=alpha<=.35:raise ValueError('Residual weight must keep a nonnegative mean component')
    return (.35-alpha)*expected+.65*poisson+alpha*empirical

def emit(x,rows,effect,own,seed,kind,prepared=None,atlas=None):
    if kind=='compact':return compact(x,rows,effect,own,seed,SETTINGS['consensus20_noise65'],prepared=prepared,atlas=atlas)
    depth,p,b=prepared[:3] if prepared is not None else moments(x)[:3];rows=np.asarray(rows,int);template=x[rows];libs=np.asarray(depth)[rows].astype(np.int64)
    np.testing.assert_array_equal(libs,np.asarray(template.sum(1)).ravel())
    factor=np.exp2(np.clip(effect,-3,3))
    if own is not None:factor[own]=.1
    wp=p*factor;wp/=wp.sum();wb=b*factor;wb/=wb.sum()
    if kind=='balanced':count,expected=balanced_counts(template,factor,seed)
    elif kind=='residual10':
        expected=libs[:,None]*wp[None,:];empirical=template.toarray().astype(float)*factor[None,:];empirical*=np.divide(libs,empirical.sum(1))[:,None]
        initial=residual_initial(expected,np.random.default_rng(seed).poisson(expected),empirical);initial/=initial.sum(1,keepdims=True)
        if atlas is None:
            from evaluate_atlas_transfer import load_atlas_model
            atlas,_=load_atlas_model()
        count=atlas.dual_moment_counts(initial,wp,wb,depths=libs,seed=int(seed),iterations=200)
    else:raise ValueError(kind)
    np.testing.assert_array_equal(count.sum(1,dtype=np.int64),libs);assert np.issubdtype(count.dtype,np.integer) and (count>=0).all()
    rp=(count/libs[:,None]).mean(0);rb=count.sum(0,dtype=np.int64)/libs.sum()
    diag={'mean_probability_L1':float(abs(rp-wp).sum()),'bulk_probability_L1':float(abs(rb-wb).sum()),'changed_transcript_fraction':float(abs(count-template.toarray()).sum()/2/libs.sum()),'global_moment_refit':kind!='balanced'}
    return sparse.csr_matrix(count),diag
