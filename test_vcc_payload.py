"""Test full-payload verification, including deliberate one-count corruption."""
import tempfile,tarfile,json
from pathlib import Path
import h5py,numpy as np,zstandard
from verify_vcc_payload import verify,compare_h5
def create(path,canonical=False):
 with h5py.File(path,'w') as f:
  g=f.create_group('X');g.attrs['shape']=np.array([2,3]);g['data']=np.array([1,3,2],np.float32 if canonical else np.int32);g['indices']=np.array([0,2,1]);g['indptr']=np.array([0,2,3])
  dt=h5py.string_dtype();o=f.create_group('obs');o.create_dataset('_index',data=np.array(['0','1'] if canonical else ['v45_0','v45_1'],object),dtype=dt);o.create_dataset('context',data=np.array(['A','B'],object),dtype=dt);o.create_dataset('target_gene',data=np.array(['X','Y'],object),dtype=dt);v=f.create_group('var');v.create_dataset('_index',data=np.array(['X','Y','Z'],object),dtype=dt)
def main():
 with tempfile.TemporaryDirectory(prefix='vcc_payload_contract_') as folder:
  p=Path(folder);src=p/'source.h5ad';dst=p/'canonical.h5ad';create(src);create(dst,True);z=p/'pred.h5ad.zst';z.write_bytes(zstandard.ZstdCompressor().compress(dst.read_bytes()));pack=p/'test.vcc'
  with tarfile.open(pack,'w') as tar:tar.add(z,arcname='pred.h5ad.zst')
  verify(pack,src,p,p/'audit.json')
  with h5py.File(dst,'r+') as f:f['X/data'][1]=4
  caught=False
  try:compare_h5(src,dst)
  except AssertionError:caught=True
  assert caught,'One corrupted count passed validation'
 from build_submission import ROOT,write_json,sha256
 write_json(ROOT/'experiments/v45_bounded_dispersion/analysis/payload_contract_test.json',{'passed':True,'checks':['complete roundtrip with official cell renumbering','integer/float32 count equality','single changed count rejected'],'script_sha256':sha256(Path(__file__))})
 print('Payload verifier contract PASS',flush=True)
if __name__=='__main__':main()

