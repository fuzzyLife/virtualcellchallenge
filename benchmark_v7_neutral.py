"""Neutral 400-cell benchmark using the full, batch-matched control population."""
import json,time
import numpy as np,pandas as pd
import build_submission as build
from fit_v7_response import raw_inputs
from evaluate_v7_response import population_statistics,matched_population,matched_rows
from evaluate_atlas_transfer import load_atlas_model
from v7_count_methods import emit
from v6_methods import moments,de
ROOT=build.ROOT
if __name__=='__main__':
    _,genes,_=build.read_spec(ROOT);x,meta,_,_=raw_inputs('H1',genes);controls=np.flatnonzero(meta.gene_target.to_numpy()=='Non-Targeting');rng=np.random.default_rng(20261008);left=[];right=[]
    for _,f in meta.iloc[controls].groupby('sample',sort=True):
        rows=rng.permutation(f.index);half=len(rows)//2;left.extend(rows[:half]);right.extend(rows[half:])
    left=np.array(left);right=np.array(right);population=population_statistics(x,meta,left);atlas,_=load_atlas_model();records=[]
    for seed in [20261008,20261019]:
        rng=np.random.default_rng(seed);reference=rng.choice(right,400,replace=False);template=matched_rows(meta,left,reference,rng);prepared=matched_population(population,meta,reference)
        actual=x[template];baseline=moments(actual)[3].sum()
        for generator in ['v4','low_noise','cal10']:
            started=time.monotonic();pred=emit(actual,np.arange(400),np.zeros(len(genes)),None,seed,generator,atlas,prepared)
            assert pred.dtype==np.int32 and (pred.data>0).all()
            if generator in ['v4','low_noise']:assert (pred!=actual).nnz==0
            result=de(pred,x[reference],genes);record={'generator':generator,'seed':seed,'false_DE_calls':int((result['q']<.05).sum()),'template_variance_ratio':float(moments(pred)[3].sum()/baseline),'seconds':time.monotonic()-started};records.append(record);print(record,flush=True)
    build.write_json(ROOT/'analysis/v7_source_inventory/neutral_population_benchmark.json',{'records':records,'matched_batches':True,'distinct_template_reference_cells':400,'population_controls':len(left),'script_sha256':build.sha256(ROOT/'benchmark_v7_neutral.py'),'scope':'Two seeds, no perturbation; full template-side population moments, 400 disjoint matched-batch reference controls.'})
