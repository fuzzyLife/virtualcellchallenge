"""Write nullable fixtures with CLI settings, read with packaging-process settings."""
import tempfile,tarfile
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad,h5py,zstandard
from scipy import sparse
from verify_v50_payload import verify,compare_h5
def main():
 original_infer=pd.options.future.infer_string
 with tempfile.TemporaryDirectory(prefix='v50_nullable_payload_test_') as folder:
  p=Path(folder);src=p/'source.h5ad';dst=p/'canonical.h5ad'
  pd.options.future.infer_string=True
  obs=pd.DataFrame({'target_gene':pd.array(['X','Y'],dtype='string'),'context':pd.Categorical(['A','B'])},index=pd.Index(pd.array(['v50_0','v50_1'],dtype='string')))
  obj=ad.AnnData(sparse.csr_matrix(np.array([[1,0,3],[0,2,0]],np.int32)),obs=obs,var=pd.DataFrame(index=['X','Y','Z']));obj.write_h5ad(src)
  obj.obs_names=pd.Index(pd.array(['0','1'],dtype='string'));obj.X=obj.X.astype(np.float32);obj.write_h5ad(dst)
  with h5py.File(dst,'r') as f:
   encoding=str(f['obs/_index'].attrs.get('encoding-type'));assert isinstance(f['obs/_index'],h5py.Group) and 'categories' not in f['obs/_index'];children=list(f['obs/_index'].keys())
  z=p/'pred.h5ad.zst';z.write_bytes(zstandard.ZstdCompressor().compress(dst.read_bytes()));pack=p/'test.vcc'
  with tarfile.open(pack,'w') as tar:tar.add(z,arcname='pred.h5ad.zst')
  pd.options.future.infer_string=False
  verify(pack,src,p,p/'audit.json')
  pd.options.future.infer_string=True
  changed=ad.read_h5ad(dst);changed.obs['target_gene']=pd.array(['Z','Y'],dtype='string');changed.write_h5ad(dst)
  pd.options.future.infer_string=False;caught=False
  try:compare_h5(src,dst)
  except AssertionError:caught=True
  assert caught,'A corrupted target label passed validation'
  pd.options.future.infer_string=True;obj.write_h5ad(dst)
  with h5py.File(dst,'r+') as f:f['X/data'][0]=2
  pd.options.future.infer_string=False;caught=False
  try:compare_h5(src,dst)
  except AssertionError:caught=True
  assert caught,'A corrupted count passed validation'
 pd.options.future.infer_string=original_infer
 import build_submission as b
 b.write_json(b.ROOT/'experiments/v50_response_prior/analysis/nullable_payload_contract_test.json',{'passed':True,'actual_index_encoding':encoding,'index_children':children,'authoring_pandas_infer_string':True,'decoding_pandas_infer_string':False,'checks':['real AnnData nullable string and categorical decoding','canonical official cell renumbering','single changed target rejected','single changed count rejected'],'script_sha256':b.sha256(Path(__file__))})
 print('V50 nullable metadata/payload tests PASS; CLI authoring=True, packaging decoder=False',flush=True)
if __name__=='__main__':main()

