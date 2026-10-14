"""Check population mean drift and integer repair; no treatment DE evaluation."""
from pathlib import Path
import json,time
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from threadpoolctl import threadpool_limits
from v22_methods import ROOT,A,raw,load_model,query_states,emit,write
from fit_v7_response import read_cache
from v16_methods import source_arrays
from v8_methods import gated
from v6_methods import moments
from v23_mean_repair import preserve_reference_means
from evaluate_atlas_transfer import load_atlas_model

def main():
 plan=json.loads((A/'plan.json').read_text());ctx=plan['contexts'][0];cell=ctx['cells'][0];data=read_cache();x,meta,mask,_=raw(ctx['family'],data['genes']);population=x[ctx['left']];_,p,bulk,_=moments(population);template=x[cell['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64);bank=source_arrays(data,[ctx['family']]);ti=int(np.flatnonzero(data['targets']==cell['target'])[0]);own=int(bank['own'][ti]);effect=gated(bank['base'][ti].copy(),p*1e6,own);model=load_model(ctx['family']);z=query_states(population,model);lookup={v:i for i,v in enumerate(ctx['left'])};z=z[[lookup[v] for v in cell['template']]];atlas,_=load_atlas_model();args=(template,effect,own,cell['seed']+ti)
 baseline,_=emit(*args,'v16',(libs,p,bulk),atlas,model,ti,z);candidate,_=emit(*args,'conditional',(libs,p,bulk),atlas,model,ti,z);t=time.monotonic();fixed,diag=preserve_reference_means(candidate.toarray(),baseline.toarray(),seed=cell['seed']+ti);diag['elapsed_seconds']=time.monotonic()-t;diag['target']=cell['target'];diag['treatment_DE_metrics_computed']=False
 out=ROOT/'experiments/v23_mean_preserving/analysis';out.mkdir(parents=True,exist_ok=True);write(out/'mean_drift_diagnostic.json',diag);print(diag,flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
