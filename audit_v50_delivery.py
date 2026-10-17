"""Final independent library/provenance alignment and package receipt audit."""
import hashlib,json
from pathlib import Path
import h5py,numpy as np,pandas as pd
pd.options.future.infer_string=False
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v50_response_prior/analysis';OUT=A.parent

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
 return h.hexdigest()

def main():
 import anndata as ad
 ready=json.loads((A/'submission_ready.json').read_text());validation=json.loads((A/'prediction_validation.json').read_text());assert ready['official_prep_passed'] and not ready['uploaded'];package=OUT/'prediction.vcc';assert package.exists() and package.stat().st_size==ready['bytes'];assert ready['package_meta']['nnz']==validation['nnz'];assert validation['n_cells']==360000 and validation['n_genes']==18533
 file=A/'prediction_control_rows.csv';frame=pd.read_csv(file);assert len(frame)==360000;np.testing.assert_array_equal(frame.original_total_counts,frame.predicted_total_counts);assert (frame.original_total_counts>0).all() and (frame.original_total_counts<=1000000).all();groups=frame.groupby(['context','target_gene'],sort=False);assert len(groups)==900 and (groups.size()==400).all() and (groups.source_control_row.nunique()==400).all();assert np.array_equal(frame.prediction_cell_offset,np.tile(np.arange(400),900))
 obj=ad.read_h5ad(OUT/'prediction.h5ad',backed='r');np.testing.assert_array_equal(obj.obs.context.to_numpy(str),frame.context.to_numpy(str));np.testing.assert_array_equal(obj.obs.target_gene.to_numpy(str),frame.target_gene.to_numpy(str));assert obj.obs_names.is_unique;obj.file.close();maximum=0;checked=0
 with h5py.File(OUT/'prediction.h5ad','r') as f:
  ptr=f['X/indptr'][:];assert len(ptr)==360001 and ptr[-1]==validation['nnz']
  for start in range(0,360000,1000):
   stop=min(start+1000,360000);lo,hi=int(ptr[start]),int(ptr[stop]);x=f['X/data'][lo:hi];totals=np.add.reduceat(x.astype(np.int64),ptr[start:stop]-lo);np.testing.assert_array_equal(totals,frame.predicted_total_counts.iloc[start:stop].to_numpy());maximum=max(maximum,int(totals.max()));checked+=len(totals)
 assert checked==360000
 receipt={'passed':True,'rows_independently_checked':checked,'groups':900,'maximum_library':maximum,'provenance_csv_sha256':sha(file),'provenance_csv_bytes':file.stat().st_size,'provenance_csv_tracked':False,'submission_ready_sha256':sha(A/'submission_ready.json'),'workspace_package_bytes':package.stat().st_size,'workspace_package_sha256_verified_by_packager':ready['sha256'],'official_archive_validation_passed':True,'checks':'Every final H5ad row library equals its recorded source-control library; obs order aligns with provenance;400 unique sampled controls per group; nonempty bounded libraries; official wrapper metadata matches validated source. Packager already verified the full remote package copy hash; not redundantly rescored.','auditor_sha256':sha(Path(__file__))}
 (A/'final_delivery_audit.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8');print('Final delivery audit passed:',checked,'rows,900 groups',flush=True)
if __name__=='__main__':main()


