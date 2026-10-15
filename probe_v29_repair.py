"""Numerical smoke on the already-spent V28 margin failure, before V29 freeze."""
from pathlib import Path
import json,time
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy import sparse
from threadpoolctl import threadpool_limits
from v26_methods import raw,fit_probabilities,integerize
from v29_pairwise_repair import repair_pairwise
from run_v26_experiment import distribution
from fit_v7_response import read_cache
from v12_cpu_metrics import de
from v22_methods import write,digest
import build_submission as b
ROOT=b.ROOT;A=ROOT/'experiments/v29_pairwise_transport/analysis'
def main():
 A.mkdir(parents=True,exist_ok=True);old=ROOT/'experiments/v28_global_transport/analysis';plan=json.loads((old/'plan.json').read_text());ctx=next(c for c in plan['contexts'] if c['family']=='CD4_Stim48_D2');case=ctx['null_cases'][-1];data=read_cache();genes=data['genes'];x,m,measured,_=raw(ctx['family'],genes);template=x[case['template']];ref=x[case['reference']];f=old/'parallel/null_CD4/CD4_Stim48_D2'/f'Non-Targeting_{case["seed"]}_v24_counts.npz';base=sparse.load_npz(f).toarray().astype(np.int32);libs=base.sum(1);p=(base/libs[:,None]).mean(0);q=fit_probabilities(template.toarray()/libs[:,None],libs,p,base.sum(0)/libs.sum());initial,n=integerize(q,libs,base.sum(0),case['seed']+270270);started=time.time();count,diag=repair_pairwise(initial,base);elapsed=time.time()-started;count=sparse.csr_matrix(count);sparse.save_npz(A/'compiled_probe_counts.npz',count);axis=np.flatnonzero(measured&~np.isin(genes,data['targets'])&(p*1e6>5));d=de(count,ref,genes);rec={'source':'Spent V28 stimulated-CD4 null margin failure, no treatment metric used','seed':case['seed'],'pairwise_code_sha256':b.sha256(ROOT/'v29_pairwise_repair.py'),'parent_plan_sha256':b.sha256(old/'plan.json'),'V24_reference_file_sha256':b.sha256(f),'count_sha256':digest(count),'seconds':elapsed,'global_repair_molecules':n,**diag,'DE_calls':int(((d['q']<.05)&d['eligible']).sum()),**distribution(count,ref,axis)};expected=json.loads((A/'numerical_probe.json').read_text());assert rec['count_sha256']==expected['count_sha256'];rec['matches_original_full_axis_digest']=True;write(A/'compiled_probe.json',rec);print(json.dumps(rec,indent=2),flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
