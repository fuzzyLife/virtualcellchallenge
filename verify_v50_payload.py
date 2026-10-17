"""V50 full CSR audit with library-backed decoding of AnnData metadata encodings."""
import json,shutil,tarfile,hashlib
from pathlib import Path
import h5py,numpy as np,zstandard,anndata as ad
def compare_h5(source,payload):
 checked={}
 with h5py.File(source,'r') as src,h5py.File(payload,'r') as dst:
  np.testing.assert_array_equal(src['X'].attrs['shape'],dst['X'].attrs['shape'])
  for key in ['data','indices','indptr']:
   a=src['X/'+key];b=dst['X/'+key];assert a.shape==b.shape
   for i in range(0,len(a),4000000):np.testing.assert_array_equal(a[i:i+4000000],b[i:i+4000000])
   checked[key]=len(a);print('Full payload matched',key,len(a),flush=True)
 left=ad.read_h5ad(source,backed='r');right=ad.read_h5ad(payload,backed='r')
 try:
  np.testing.assert_array_equal(right.obs_names.to_numpy(str),np.arange(left.n_obs).astype(str))
  for key in ['context','target_gene']:np.testing.assert_array_equal(left.obs[key].to_numpy(str),right.obs[key].to_numpy(str))
  np.testing.assert_array_equal(left.var_names.to_numpy(str),right.var_names.to_numpy(str));checked['axis_fields']=4
 finally:left.file.close();right.file.close()
 return checked
def verify(package,staged,scratch,receipt_path):
 payload=Path(scratch)/'roundtrip_payload.h5ad';assert payload.parent==Path(scratch)
 with tarfile.open(package,'r:*') as tar:
  matches=[m for m in tar.getmembers() if m.name=='pred.h5ad.zst'];assert len(matches)==1 and matches[0].isfile();stream=tar.extractfile(matches[0]);assert stream is not None
  with zstandard.ZstdDecompressor().stream_reader(stream) as reader,payload.open('wb') as out:shutil.copyfileobj(reader,out,8*1024*1024)
 result={'passed':True,'member':'pred.h5ad.zst','full_CSR_comparison':compare_h5(staged,payload),'metadata_decoder':'anndata.read_h5ad(backed=r), supporting datasets, categorical and nullable strings','auditor_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'All CSR values/indices/pointers and gene/context/target labels exactly match the staged source. Official sequential cell renumbering checked explicitly. Wrapper-only validation is supplemented by this full payload comparison.'}
 Path(receipt_path).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');payload.unlink();return result

