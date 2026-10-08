"""Test upstream-inspired centering on held-out CD4 targets and dual moments on control reconstruction."""
from pathlib import Path
import argparse
import importlib.util
import json
import sys
import time
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as build
import fit_lineage_submission as v4
import fit_context_ensemble as v3
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'analysis/upstream_review'
FIELDS=['loss','zero_loss','selection_score','direction_fidelity_proxy','jaccard_proxy','mean_predicted_calls','mean_reference_calls']


def source_center(values,valid,conditions,train):
    centers={}
    for condition in v4.CONDITIONS:
        use=train&(conditions==condition);mask=valid[use];n=mask.sum(axis=0)
        centers[condition]=np.divide(np.where(mask,values[use],0).sum(axis=0),n,out=np.zeros(values.shape[1]),where=n>0)
    return np.stack([centers[c] for c in conditions])


def centered_prediction(values,errors,center,setting):
    return np.clip(setting['gain']*(values-setting['centering']*center)*.04/(.04+errors**2),-3,3)


def aggregate(records):
    n=sum(x['contrasts'] for x in records)
    result={field:sum(x[field]*x['contrasts'] for x in records)/n for field in FIELDS}
    result.update({'contrasts':n,'relative_zero_loss':result['loss']/result['zero_loss']})
    return result


def select(records):
    baseline=next(x for x in records if x['setting']=={'gain':1.,'centering':0.})
    # For this error-focused follow-up, training error may not exceed the V4 fixed baseline.
    valid=[x for x in records if x['relative_zero_loss']<=min(1.02,baseline['relative_zero_loss'])+1e-12]
    return max(valid,key=lambda x:(x['selection_score'],-x['setting']['centering'],-x['setting']['gain']))


def centering():
    rows,left,right,lq,rq,valid,le,re,conditions,targets,sl,sr,names=v4.donor_inputs()
    folds=v3.target_folds(targets);grid=[{'gain':gain,'centering':center} for gain in [.5,.75,1.] for center in [0.,.5,1.]]
    outer=[];nested=[]
    def predictions(train):
        return source_center(left,valid,conditions,train),source_center(right,valid,conditions,train)
    for fold in range(5):
        test=folds==fold;train=~test;inner=v3.target_folds(targets[train],3);indices=np.flatnonzero(train)
        scores=[[] for _ in grid]
        for part in range(3):
            tune=np.zeros(len(rows),bool);tune[indices[inner==part]]=True;fit=train&~tune
            lc,rc=predictions(fit)
            for ci,setting in enumerate(grid):
                pl=centered_prediction(left,le,lc,setting);pr=centered_prediction(right,re,rc,setting)
                scores[ci].append(v4.evaluate(pl,pr,left,right,lq,rq,valid,sl,sr,tune))
        choices=[{'candidate':i,'setting':setting,**aggregate(scores[i])} for i,setting in enumerate(grid)]
        chosen=select(choices);lc,rc=predictions(train)
        for ci,setting in enumerate(grid):
            pl=centered_prediction(left,le,lc,setting);pr=centered_prediction(right,re,rc,setting)
            result=v4.evaluate(pl,pr,left,right,lq,rq,valid,sl,sr,test)
            record={'fold':fold,'candidate':ci,'setting':setting,**result};outer.append(record)
            if ci==chosen['candidate']:nested.append({**record,'train_targets':sorted(set(targets[train])),'test_targets':sorted(set(targets[test]))})
        print('Centering fold',fold,chosen['setting'],'held-out relative error',nested[-1]['relative_zero_loss'],flush=True)
    candidates=[{'candidate':i,'setting':setting,**aggregate([r for r in outer if r['candidate']==i])} for i,setting in enumerate(grid)]
    baseline=next(r for r in candidates if r['setting']=={'gain':1.,'centering':0.});procedure=aggregate(nested)
    accepted=procedure['selection_score']>baseline['selection_score']+1e-6 and procedure['loss']<=baseline['loss'] and any(r['setting']['centering']>0 for r in nested)
    for label,records in [('centering_candidates',candidates),('centering_outer_scores',outer)]:
        pd.DataFrame([{**x,'setting':json.dumps(x['setting'],sort_keys=True)} for x in records]).to_csv(OUT/(label+'.csv'),index=False)
    build.write_json(OUT/'centering_validation.json',{'v4_fixed_baseline':baseline,'nested_procedure':procedure,'nested_folds':nested,
        'centering_change_supported':bool(accepted),'targets':len(set(targets)),'contrasts':len(rows),'genes':left.shape[1],
        'grid':grid,'source_donor_groups':names,'script_sha256':build.sha256(Path(__file__)),
        'policy':'Five outer target folds and three inner target folds; all target conditions kept together; own-gene-invalid coordinates excluded from centers. Candidate training error cannot exceed fixed V4; acceptance also requires no held-out loss increase and improved composite proxy.',
        'scope':'Unweighted source centering inspired by AtlasShift, tested on cached CD4 donor contrasts; raw expression-error score and anonymous-context transfer are not observed here. Centering may remove genuine shared within-lineage biology.'})
    print('Centering supported:',accepted,'baseline',baseline['selection_score'],'nested',procedure['selection_score'],flush=True)


def load_atlas_model():
    path=ROOT/'Virtual-Cell-Challenge-2026/model.py';review=json.loads((OUT/'repository_review.json').read_text())
    repo=next(r for r in review['repositories'] if r['repository']=='Virtual-Cell-Challenge-2026')
    expected=next(r['sha256'] for r in repo['tracked_files'] if r['path']=='model.py')
    assert build.sha256(path)==expected
    spec=importlib.util.spec_from_file_location('reviewed_atlas_model',path);module=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=module;spec.loader.exec_module(module)
    return module,expected


def probability_stats(x):
    library=np.asarray(x.sum(axis=1)).ravel();assert (library>0).all()
    normalized=x.multiply((1/library)[:,None]).tocsr()
    mean=np.asarray(normalized.mean(axis=0)).ravel();bulk=np.asarray(x.sum(axis=0)).ravel();bulk=bulk/bulk.sum()
    variance=np.maximum(np.asarray(normalized.power(2).mean(axis=0)).ravel()-mean**2,0)
    return library,mean,bulk,variance


def moments():
    model,model_sha=load_atlas_model();records=[]
    for ci,context in enumerate('ABC'):
        original=ad.read_h5ad(ROOT/f'context_{context}.h5ad');x=original.X.tocsr().astype(float)
        guides=original.obs.ntc_id.astype(str).to_numpy();rng=np.random.default_rng(20261005+ci)
        order=rng.permutation(len(guides));training=order[:len(order)//2];test=order[len(order)//2:]
        libs,p,bulk,_=probability_stats(x[training]);_,test_p,test_bulk,test_var=probability_stats(x[test])
        donor=build.guide_stratified_sample(guides[training],libs,400,rng);baseline=x[training[donor]].toarray().astype(np.int32)
        def describe(counts,label,pool):
            counts=sparse.csr_matrix(counts);depth,cm,bm,var=probability_stats(counts)
            return {'context':context,'method':label,'pool_k':pool,'cells':counts.shape[0],
                'mean_CPM_probability_L1_from_training':float(np.abs(cm-p).sum()),
                'pooled_probability_L1_from_training':float(np.abs(bm-bulk).sum()),
                'held_out_control_log_bulk_MSE':float(np.mean((np.log1p(50000*bm)-np.log1p(50000*test_bulk))**2)),
                'held_out_control_log_mean_CPM_MSE':float(np.mean((np.log1p(50000*cm)-np.log1p(50000*test_p))**2)),
                'CPM_variance_ratio_to_heldout_controls':float(var.sum()/test_var.sum()),
                'mean_detected_genes':float(counts.getnnz(axis=1).mean()),'mean_library':float(depth.mean()),
                'library_CV':float(depth.std()/depth.mean()),'positive_integer_counts':bool(np.issubdtype(counts.dtype,np.integer) and (counts.data>0).all())}
        records.append(describe(baseline,'original_control_sample',1))
        for pool in [1,4]:
            selected=donor.copy() if pool==1 else rng.choice(len(training),400*pool,replace=False)
            selected=selected[np.argsort(libs[selected],kind='stable')]
            template=x[training[selected]].toarray()/libs[selected,None]
            template=template.reshape(400,pool,x.shape[1]).mean(axis=1)
            depths=np.rint(libs[selected].reshape(400,pool).mean(axis=1)).astype(np.int64)
            started=time.monotonic()
            with threadpool_limits(limits=2):
                counts=model.dual_moment_counts(template,p,bulk,depths=depths,seed=20261005+ci*10+pool)
            np.testing.assert_array_equal(counts.sum(axis=1),depths)
            record=describe(counts,'dual_moment_control_reconstruction',pool);record['seconds']=time.monotonic()-started;records.append(record)
            print(context,'dual moments, pool',pool,'L1',record['pooled_probability_L1_from_training'],'variance ratio',record['CPM_variance_ratio_to_heldout_controls'],flush=True)
        original.file.close();del original,x
    pd.DataFrame(records).to_csv(OUT/'control_moment_diagnostics.csv',index=False)
    build.write_json(OUT/'control_moment_scope.json',{'upstream_model_sha256':model_sha,'script_sha256':build.sha256(Path(__file__)),
        'training_test_cells_per_context':9200,'predicted_cells_each_group':400,'upstream_pool_sizes':[1,4],'pool_one_uses_same_control_donors_as_baseline':True,
        'all_groups_integer_and_preserve_cell_depth':True,
        'scope':'No-response control reconstruction, distinct training/test controls; no perturbation labels, own-target knockdown or hidden response truth. Tests moment geometry and dispersion, not leaderboard accuracy. Do not automatically replace observed-cell variability with pooled cells.'})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--centering-only',action='store_true');parser.add_argument('--moments-only',action='store_true');args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    if not args.moments_only:centering()
    if not args.centering_only:moments()
