"""Post-freeze descriptive localization inside V27 anchor projection.

These are numerical diagnostics on the spent control case, not additional
eligible candidates or a revision of the locked V27 release rules.
"""
from pathlib import Path
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy import sparse
from threadpoolctl import threadpool_limits
from v27_methods import ROOT,A,load_plan,write,digest
from v26_methods import raw,fit_probabilities,integerize,preserve_depth_margins,repair_cpm,moments
from fit_v7_response import read_cache
from v12_cpu_metrics import de
from run_v26_experiment import distribution
import build_submission as b

def main():
 plan=load_plan();data=read_cache();genes=data['genes'];case=plan['contexts'][0]['null_cases'][0];x,m,measured,_=raw('H1',genes);template=x[case['template']];ref=x[case['reference']];base=sparse.load_npz(A/'sampling/v24_reference.npz').toarray().astype(np.int32);libs=base.sum(1);p=(base/libs[:,None]).mean(0);bulk=base.sum(0)/libs.sum();initial=template.toarray().astype(float)/libs[:,None];prob=fit_probabilities(initial,libs,p,bulk);seed=case['seed']+270270;global_count,n=integerize(prob,libs,base.sum(0),seed);axis=np.flatnonzero(measured&~np.isin(genes,data['targets'])&(p*1e6>5));records=[];folder=A/'projection_diagnostic';folder.mkdir(exist_ok=True)
 def record(name,v):
  c=sparse.csr_matrix(v);d=de(c,ref,genes);sparse.save_npz(folder/(name+'.npz'),c);row={'stage':name,'DE_calls':int(((d['q']<.05)&d['eligible']).sum()),'mean_reference_L1':float(abs((v/libs[:,None]).mean(0)-p).sum()),'count_sha256':digest(c),**distribution(c,ref,axis)};records.append(row);print('V27 anchor projection',row,flush=True)
 record('global_integer_repair',global_count)
 no_depth,k=repair_cpm(global_count.copy(),base);record('global_CPM_without_strata',no_depth)
 with_depth,nn=preserve_depth_margins(global_count,base,seed+1);record('depth_strata_repair',with_depth)
 final,kk=repair_cpm(with_depth,base);record('frozen_final_CPM',final);assert digest(sparse.csr_matrix(final))==digest(sparse.load_npz(A/'sampling/raw_anchor_projected.npz'))
 pd.DataFrame(records).to_csv(A/'projection_diagnostic.csv',index=False);write(A/'projection_diagnostic.json',{'plan_sha256':b.sha256(A/'plan.json'),'scores_sha256':b.sha256(A/'projection_diagnostic.csv'),'diagnostician_sha256':b.sha256(Path(__file__)),'global_repair_molecules':n,'no_depth_CPM_molecules':k,'depth_molecules':nn,'final_CPM_molecules':kk,'reproduces_frozen_anchor_digest':True,'scope':'Post-freeze diagnostic of spent null only. No additional eligible arm or treatment inference; no model, criterion or decision revision.'})
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
