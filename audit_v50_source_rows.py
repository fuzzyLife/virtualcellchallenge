"""Independently tie every prediction provenance row back to the supplied controls."""
from pathlib import Path
import json,hashlib
import h5py,numpy as np,pandas as pd,anndata as ad
import build_submission as b
ROOT=b.ROOT;A=ROOT/'experiments/v50_response_prior/analysis'
def main():
 frame=pd.read_csv(A/'prediction_control_rows.csv',dtype={'source_control_cell_id':str,'ntc_id':str,'context':str,'target_gene':str});assert len(frame)==360000;records=[]
 manifest,genes,targets=b.read_spec(ROOT)
 for context in manifest['contexts']:
  path=ROOT/f'context_{context}.h5ad';obj=ad.read_h5ad(path,backed='r');ids=obj.obs_names.to_numpy(str);guides=obj.obs.ntc_id.to_numpy(str);np.testing.assert_array_equal(obj.var_names.to_numpy(str),genes);obj.file.close()
  with h5py.File(path,'r') as f:
   assert f['X'].attrs['encoding-type']=='csr_matrix';ptr=f['X/indptr'][:];assert (np.diff(ptr)>0).all();libs=np.zeros(len(ptr)-1,np.int64)
   for start in range(0,len(libs),500):
    stop=min(start+500,len(libs));lo=int(ptr[start]);hi=int(ptr[stop]);values=f['X/data'][lo:hi];assert np.isfinite(values).all() and (values>=0).all() and np.equal(values,np.floor(values)).all()
    libs[start:stop]=np.add.reduceat(values.astype(np.int64),ptr[start:stop]-lo)
  sub=frame[frame.context==context];rows=sub.source_control_row.to_numpy(np.int64);assert (rows>=0).all() and (rows<len(libs)).all()
  np.testing.assert_array_equal(ids[rows],sub.source_control_cell_id.to_numpy(str));np.testing.assert_array_equal(guides[rows],sub.ntc_id.to_numpy(str));np.testing.assert_array_equal(libs[rows],sub.original_total_counts.to_numpy(np.int64));np.testing.assert_array_equal(libs[rows],sub.predicted_total_counts.to_numpy(np.int64))
  records.append({'context':context,'prediction_rows_verified':len(sub),'source_rows':len(libs),'source_sha256':b.sha256(path),'source_bytes':path.stat().st_size});print('Original source rows verified',context,len(sub),flush=True)
 b.write_json(A/'source_row_audit.json',{'passed':True,'provenance_csv_sha256':b.sha256(A/'prediction_control_rows.csv'),'rows':sum(r['prediction_rows_verified'] for r in records),'records':records,'checks':'Independent streamed original H5ad libraries, source barcodes and NTC identifiers match every selected prediction/provenance row. Gene order matches declared competition axis.','script_sha256':b.sha256(Path(__file__))})
if __name__=='__main__':main()


