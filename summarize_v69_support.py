"""Describe how much source-trained dispersion evidence reaches each public case."""
import numpy as np,pandas as pd
from run_v35_v36 import contexts,key
from run_v69_experiment import unique
from v62_conserved_response import data
import v69_count_bridge as m
from v57_pie import read,write,sha


def main():
    d=data();rows=[];banks={}
    for ctx in contexts():
        group=key(ctx['family'])
        if group not in banks:banks[group]=m.pooled((group,))[0]
        p=m.ROOT/f"experiments/v66_assay_strength/analysis/predictions/{ctx['family']}.npz"
        with np.load(p,allow_pickle=False) as z:usable=z['measured']&(np.nan_to_num(z['query'])>0)
        for i in unique(ctx):
            target=ctx['cells'][i]['target'];ti=list(d['targets']).index(target);delta=banks[group][ti];supported=(delta!=0)&usable
            rows.append(dict(context=ctx['family'],study=group,target=target,case_index=i,source_supported_genes=int(np.count_nonzero(delta)),query_usable_supported_genes=int(supported.sum()),mean_absolute_transferred_log_alpha=float(np.mean(abs(.5*delta[usable]))),any_source_shape_evidence=bool(supported.any())))
    p=m.A/'shape_support_by_case.csv';pd.DataFrame(rows).to_csv(p,index=False)
    write(m.A/'shape_support_summary.json',dict(cases=len(rows),cases_with_any_source_shape_evidence=sum(r['any_source_shape_evidence'] for r in rows),csv_sha256=sha(p),scope='Descriptive source coverage, not a performance metric. Source public studies excluded; no recipient treatments read.'))
    print(read(m.A/'shape_support_summary.json'),flush=True)

if __name__=='__main__':main()
