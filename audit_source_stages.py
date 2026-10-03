"""Measure source suppression and matched-size non-targeting pseudo-perturbations.

Null controls diagnose source-noise specificity. They do not validate biological
cross-context transfer, and no candidate is selected by hidden leaderboard labels.
"""
from pathlib import Path
import argparse
import json
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import build_submission as build


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--null-groups',type=int,default=40)
    args=parser.parse_args()
    root=args.root.resolve(); output=root/'analysis/model_audit'; output.mkdir(exist_ok=True)
    _,genes,targets=build.read_spec(root)
    with h5py.File(root/'reference/K562_gwps_raw_bulk_01.h5ad','r') as f:
        labels=f['obs/gene_transcript'].asstr()[:]; cells=f['obs/num_cells_filtered'][:].astype(float)
        fold=f['obs/fold_expr'][:]; core=f['obs/core_control'][:].astype(bool)
        field=f['var/gene_name']; symbols=f[field.attrs['categories']].asstr()[:][field[:]]
        matrix=f['X'][:].astype(float)
    source_targets=np.array([s.split('_',2)[1] if s.split('_',1)[0].isdigit() else '' for s in labels])
    valid=np.isfinite(cells)&(cells>=30); control=core&(source_targets=='non-targeting')&valid
    shared=sorted(set(symbols)&set(genes)); lookup={g:i for i,g in enumerate(shared)}
    mask=np.array([g in lookup for g in symbols])
    mapping=sparse.csr_matrix((np.ones(mask.sum()),(np.flatnonzero(mask),[lookup[g] for g in symbols[mask]])),shape=(len(symbols),len(shared)))
    matrix=mapping.T.dot(matrix.T).T
    basal=np.average(matrix[control],axis=0,weights=cells[control])
    source_control_indices=np.flatnonzero(control)
    source_control_n=cells[control].sum()
    groups=[]
    for target in targets:
        rows=np.flatnonzero(valid&(source_targets==target)&np.isfinite(fold)&(fold<=.5))
        if len(rows):groups.append((target,rows))
    rows_out=[]; n_by_target=[]
    def evaluate(label, indices, base, base_n, kind):
        mean=np.average(matrix[indices],axis=0,weights=cells[indices]); n=cells[indices].sum()
        cpm=mean*1e6/mean.sum(); base_cpm=base*1e6/base.sum()
        raw=np.log2((cpm+1)/(base_cpm+1))
        noise=np.sqrt((1/np.maximum(n*mean,1)+1/np.maximum(base_n*base,1))/(np.log(2)**2))
        downstream=np.array(shared)!=label
        estimates={}
        for name,z,reliability in [('raw',0,False),('reliability_only',0,True),('original_threshold',1.96,True),('less_threshold',.5,True)]:
            estimate=np.sign(raw)*np.maximum(np.abs(raw)-z*noise,0)
            if reliability:estimate=np.clip(estimate,-3,3)*n/(n+50)
            estimates[name]=estimate
            rows_out.append({'kind':kind,'stage':name,'target_gene':label,'source_cells':int(n),
                             'nonzero_downstream_genes':int(np.count_nonzero(estimate[downstream])),
                             'genes_gt_point1_at_shrinkage_point3':int((np.abs(.3*estimate[downstream])>.1).sum()),
                             'genes_gt_point1_at_full_strength':int((np.abs(estimate[downstream])>.1).sum())})
        return estimates
    cached=np.load(root/'analysis/reference_model_arrays.npz')
    cached_effects=cached['effects'][0]
    cached_targets={target:i for i,target in enumerate(targets)}
    maximum_difference=0.
    for target,indices in groups:
        estimates=evaluate(target,indices,basal,source_control_n,'perturbed')
        original=cached_effects[cached_targets[target]]
        maximum_difference=max(maximum_difference,float(np.abs(estimates['original_threshold']-original).max()))
        n_by_target.append(cells[indices].sum())
    # Independent control groups versus the remaining controls, with approximately
    # the same cell counts as source perturbations. Entire guide profiles are used.
    rng=np.random.default_rng(20261003)
    selected_n=rng.choice(n_by_target,size=args.null_groups,replace=True)
    for i,desired_n in enumerate(selected_n):
        permutation=rng.permutation(source_control_indices)
        cumulative=cells[permutation].cumsum()
        take=int(np.argmin(np.abs(cumulative-desired_n)))+1
        indices=permutation[:take]
        remaining=permutation[take:]
        baseline=np.average(matrix[remaining],axis=0,weights=cells[remaining])
        evaluate(f'null-{i}',indices,baseline,cells[remaining].sum(),'non_targeting_null')
    assert maximum_difference<1e-6, 'Source reconstruction differs from submitted reference cache'
    frame=pd.DataFrame(rows_out)
    frame.to_csv(output/'source_stage_and_nulls.csv',index=False)
    summary=frame.groupby(['kind','stage'])[['source_cells','nonzero_downstream_genes','genes_gt_point1_at_shrinkage_point3','genes_gt_point1_at_full_strength']].agg(['median','mean','min','max'])
    summary.columns=['_'.join(column) for column in summary.columns]
    summary.reset_index().to_csv(output/'source_stage_and_nulls_summary.csv',index=False)
    build.write_json(output/'source_stage_checks.json',{'seed':20261003,'null_groups':args.null_groups,
                     'covered_perturbations':len(groups),'maximum_original_cache_difference':maximum_difference,
                     'null_design':'Disjoint groups of NTC guide pseudobulk profiles versus remaining NTC guides; approximate perturbation cell-count matching.',
                     'limitation':'Source specificity diagnostic only; different constructs can share systematic batch effects. No cross-context validation.'})
    print(summary.to_string(),flush=True)


if __name__=='__main__':main()
