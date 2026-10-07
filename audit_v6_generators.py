"""Paired actual-context V4/V5 response/generator ablations; activity, not accuracy."""
import json
import argparse
import gc
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from audit_prediction import normalized_counts,ControlRankTest
from evaluate_atlas_transfer import load_atlas_model
from v6_methods import moments,emit,CANDIDATES
from build_v5_submission import desired_profiles
import build_submission as build
ROOT=build.ROOT;A=ROOT/'experiments/v6_fidelity_calibrated/analysis'


def main(partial_only=False):
    _,genes,targets=build.read_spec(ROOT);models=[]
    for folder in ['v4_lineage_response','v5_atlas_calibrated']:
        with np.load(ROOT/f'experiments/{folder}/analysis/reference_model_arrays.npz',allow_pickle=False) as z:
            models.append({k:z[k].copy() for k in z.files})
    atlas,_=load_atlas_model();names=['AKT2','SMARCA5','VIM','S100A11','FXYD5','KIF21B'];records=[]
    for ci,context in enumerate('ABC'):
        manifest,_,_=build.read_spec(ROOT);x,guides,library,_,_,_=build.load_control(ROOT,context,genes,manifest)
        _,p,b,var=moments(x);prepared=(library,p,b,np.asarray(x.mean(axis=0)).ravel())
        norm=normalized_counts(x);mean=np.asarray(norm.mean(axis=0)).ravel();columns=np.flatnonzero(mean>5);rank=ControlRankTest(norm[:,columns])
        for target in names:
            ti=int(np.flatnonzero(targets==target)[0]);own=int(np.flatnonzero(genes==target)[0]);seed=20261007+ci*100000+ti
            rows=build.guide_stratified_sample(guides,library,400,np.random.default_rng(np.random.SeedSequence([20261007,ci,ti])))
            for name in (['v4_cal10','v4_cal25'] if partial_only else ['v4','v4_dual','v5_response_v4_generator','v5_dual']):
                effect=models[1 if name.startswith('v5') else 0]['effects'][ci,ti]
                requested=desired_profiles(p,b,effect,models[1]['bulk_delta'][ci,ti],own)[1] if name=='v5_dual' else None
                generated,diagnostic=emit(x,rows,effect,own,seed,CANDIDATES[name],atlas,requested,prepared)
                pv,q,lfc=rank.test(normalized_counts(generated)[:,columns]);called=(q<.05)&(genes[columns]!=target)
                _,mp,mb,gv=moments(generated)
                records.append({'context':context,'target_gene':target,'candidate':name,'predicted_DE_calls':int(called.sum()),
                    'positive_calls':int((called&(lfc>0)).sum()),'negative_calls':int((called&(lfc<0)).sum()),
                    'variance_ratio_to_controls':float(gv.sum()/var.sum()),'mean_detected_genes':float(generated.getnnz(axis=1).mean()),**diagnostic})
            print(context,target,'paired 2x2 generator audit complete',flush=True)
        del x,norm,rank;gc.collect()
    frame=pd.DataFrame(records)
    if partial_only:
        previous=pd.read_csv(A/'context_generator_ablation.csv')
        previous=previous[~previous.candidate.isin(['v4_cal10','v4_cal25'])]
        frame=pd.concat([previous,frame],ignore_index=True)
    frame.to_csv(A/'context_generator_ablation.csv',index=False)
    build.write_json(A/'context_ablation_scope.json',{'groups':len(frame),'targets':names,'paired_original_control_draws':True,
        'scope':'Generated response activity against provided controls. No hidden perturbation truth, hence no direction-fidelity accuracy attribution.',
        'methods_sha256':build.sha256(ROOT/'v6_methods.py'),'script_sha256':build.sha256(__import__('pathlib').Path(__file__))})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--partial-only',action='store_true');args=parser.parse_args()
    with threadpool_limits(limits=2):main(args.partial_only)
