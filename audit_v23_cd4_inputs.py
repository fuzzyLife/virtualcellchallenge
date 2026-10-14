"""Independently audit raw-CD4 identities/lanes and sampled aligned counts."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
pd.options.future.infer_string=False
from scipy import sparse
from v22_methods import ROOT,A
import build_submission as b

def main():
 folder=A/'cd4_selected';rec=json.loads((folder/'extraction.json').read_text());d=pd.read_csv(folder/'cells.csv');axis=pd.read_csv(folder/'genes.csv');x=sparse.load_npz(folder/'counts.npz');source=Path(rec['source']);obj=ad.read_h5ad(source,backed='r');assert (d.low_quality==False).all();np.testing.assert_array_equal(d['sample'].astype(str),d.lane_id.astype(str));np.testing.assert_array_equal(d.cell_id,obj.obs_names.to_numpy()[d.source_row.to_numpy(int)]);names=obj.obs.perturbation.astype(str).to_numpy();expected=names[d.source_row.to_numpy(int)].copy();expected[np.isin(expected,['control'])]='Non-Targeting';np.testing.assert_array_equal(expected,d.gene_target.to_numpy(str));np.testing.assert_array_equal(np.asarray(x.sum(1)).ravel(),d.total_counts.to_numpy());assert np.issubdtype(x.dtype,np.integer) and (x.data>=0).all();ids=np.sort(np.random.default_rng(23093).choice(len(d),128,replace=False));rows=d.iloc[ids].source_row.to_numpy(int);original=obj.X[rows].tocsr();symbols=obj.var_names.to_numpy(str);lookup={g:i for i,g in enumerate(axis.gene_name.to_numpy(str))};dest=np.array([lookup.get(g,-1) for g in symbols]);coo=original.tocoo();keep=dest[coo.col]>=0;aligned=sparse.csr_matrix((coo.data[keep],(coo.row[keep],dest[coo.col[keep]])),shape=(len(ids),len(axis)));assert (aligned!=x[ids]).nnz==0;obj.file.close();reserve=json.loads((ROOT/'experiments/v23_mean_preserving/analysis/reserve_plan.json').read_text());assert not set(reserve['left'])&set(reserve['right']);assert (d.gene_target.iloc[reserve['left']]=='Non-Targeting').all() and (d.gene_target.iloc[reserve['right']]=='Non-Targeting').all();assert len(reserve['cells'])==8
 for cell in reserve['cells']:
  assert len(set(cell['actual']))==len(set(cell['template']))==400;assert (d.gene_target.iloc[cell['actual']]==cell['target']).all();assert set(cell['template'])<=set(reserve['left']);assert not set(cell['actual'])&(set(reserve['left'])|set(reserve['right']))
 b.write_json(A/'cd4_raw_independent_audit.json',{'passed':True,'all_retained_cell_ids_labels_libraries_and_lane_names_verified':len(d),'all_selected_cells_pass_original_quality_flag':True,'random_original_rows_with_every_gene_count_rechecked':128,'every_sampled_aligned_count_matches_original':True,'reserve_real_cells_verified':3200,'disjoint_control_left':len(reserve['left']),'disjoint_control_right':len(reserve['right']),'source_md5':rec['source_md5'],'extraction_sha256':b.sha256(folder/'extraction.json'),'auditor_sha256':b.sha256(Path(__file__))});print('CD4 raw identities/lanes and128 complete count rows independently verified',flush=True)
if __name__=='__main__':main()
