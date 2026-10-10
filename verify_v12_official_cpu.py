"""Run the pinned cell-eval2 CPU backend independently in Python 3.13."""
import json,platform,importlib.metadata as md,hashlib
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from cell_eval2.de_compute import compute_de
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v12_distribution_calibrated/analysis'
def main():
 with np.load(A/'cpu_parity_fixture.npz') as z:genes=z['genes'];c=z['control'];t=z['treated']
 data=ad.AnnData(sparse.csr_matrix(np.vstack([c,t])),obs=pd.DataFrame({'target':['non-targeting']*len(c)+['TARGET']*len(t)},index=[str(i) for i in range(len(c)+len(t))]),var=pd.DataFrame(index=genes))
 frame=compute_de(data,backend='scanpy',groupby='target',reference='non-targeting',mean_calc='arithmetic',epsilon=1e-9,input_type='counts',target_sum=1e6,clip_value=None,filter_gene_min_cpm_cell=5.,fdr_scope='per_pert',threads=2,device='cpu')
 frame.write_csv(A/'official_cpu_de.csv');f=pd.DataFrame(frame.to_dict(as_series=False)).set_index('feature').reindex(genes)
 with np.load(A/'cpu_parity_local.npz') as z:
  eligible=z['eligible'];assert set(frame['feature'])==set(genes[eligible]);differences={}
  for key,column in [('p','p_value'),('q','p_adj'),('lfc','log2_fold_change')]:
   first=z[key][eligible];second=f[column].to_numpy()[eligible]
   np.testing.assert_allclose(first,second,rtol=1e-6,atol=1e-9,err_msg=key)
   differences[key]=float(np.max(abs(first-second)))
 versions={d:md.version(d) for d in ['cell-eval2','scanpy','numpy','scipy','anndata','polars','pandas']}
 result={'passed':True,'backend':'scanpy CPU','official_commit':'0ce25cf4c495e76e60dfca4ddc8f822507aa0e18','python':platform.python_version(),'versions':versions,'eligible_genes':int(eligible.sum()),'max_absolute_difference':differences,
 'scope':'Independent full official actual-count DE path. Metric-layer fixtures tested separately. Does not assert hosted GPU backend parity.'}
 (A/'official_cpu_parity.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result),flush=True)
if __name__=='__main__':main()
