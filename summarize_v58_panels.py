"""Descriptive common-control panel scoring of already frozen V57 count tests."""
import numpy as np
import pandas as pd
from v58_pie import *
from fit_v58_pie_decoder import inputs
from v25_methods import raw
from audit_transfer_validation import pds,profile_jackknife_chunked,corrected_expression

def main():
 pp=plan();rows=[]
 for source in SOURCES:
  z=inputs(source);genes=z['genes'];x,meta,measured,receipt=raw(source,genes);ntc=np.flatnonzero(meta.gene_target.to_numpy(str)=='Non-Targeting');cp,cj=profile_jackknife_chunked(x[ntc]);axis=~np.isin(genes,z['targets'])
  seeds=pp['config']['count_seeds'] if source.startswith('CD4') else pp['config']['count_seeds'][:1]
  for seed in seeds:
   samples=[];targets=pp['cohorts'][source]['count_pilot']
   for target in targets:
    file=CACHE/'count_statistics'/f'{source}_{target}_{seed}.npz';rec=read(A/'counts'/file.with_suffix('.json').name);assert sha(file)==rec['statistics_sha256']
    with np.load(file,allow_pickle=False) as v:samples.append({k:v[k].copy() for k in v.files})
   truth=np.stack([s['truth_profile'] for s in samples]);tj=np.stack([s['truth_jackknife'] for s in samples]);own=np.array([list(genes).index(t) for t in targets])
   np.testing.assert_allclose(pds(np.broadcast_to(cp,truth.shape)[:,axis],truth[:,axis],cp[axis]),.5,atol=1e-12)
   for arm in pp['config']['arms']:
    pred=np.stack([s[arm+'_profile'] for s in samples]);pj=np.stack([s[arm+'_jackknife'] for s in samples]);ex=corrected_expression(pred,truth,cp,pj,tj,cj,own)
    rows.append({'source':source,'seed':seed,'targets':len(targets),'arm':arm,'pds':pds(pred[:,axis],truth[:,axis],cp[axis]),**ex})
  del x
 frame=pd.DataFrame(rows);frame.to_csv(A/'panel_metrics.csv',index=False)
 write(A/'panel_metrics_receipt.json',{'passed':True,'script_sha256':sha(Path(__file__)),'plan_sha256':sha(A/'plan.json'),'sha256':sha(A/'panel_metrics.csv'),'scope':'Secondary descriptive common-control panel calculation. One pooled NTC profile per context, same for every target and arm. Four/eight target micro-panels do not establish 300-target PDS. No selection or invented overall score. Expression numerators and denominators retained; unstable ratios not averaged into a ranking forecast.'})
 print(frame.to_string(index=False),flush=True)

if __name__=='__main__':main()
