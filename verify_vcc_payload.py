"""Stream one official archive member and compare ALL CSR values and axes."""
import json,shutil,tarfile
from pathlib import Path
import h5py,numpy as np,zstandard
def compare_h5(source,payload):
 checked={}
 with h5py.File(source,'r') as src,h5py.File(payload,'r') as dst:
  np.testing.assert_array_equal(src['X'].attrs['shape'],dst['X'].attrs['shape'])
  for key in ['data','indices','indptr']:
   a=src['X/'+key];b=dst['X/'+key];assert a.shape==b.shape
   for i in range(0,len(a),4000000):np.testing.assert_array_equal(a[i:i+4000000],b[i:i+4000000])
   checked[key]=len(a);print('Full payload matched',key,len(a),flush=True)
  def labels(f,path):
   x=f[path]
   if isinstance(x,h5py.Group):return x['categories'].asstr()[:][x['codes'][:]]
   return x.asstr()[:]
  np.testing.assert_array_equal(labels(dst,'obs/_index'),np.arange(src['X'].attrs['shape'][0]).astype(str))
  for path in ['obs/context','obs/target_gene','var/_index']:np.testing.assert_array_equal(labels(src,path),labels(dst,path))
  checked['axis_fields']=4
 return checked
def verify(package,staged,scratch,receipt_path):
 payload=Path(scratch)/'roundtrip_payload.h5ad';assert payload.parent==Path(scratch)
 with tarfile.open(package,'r:*') as tar:
  matches=[m for m in tar.getmembers() if m.name=='pred.h5ad.zst'];assert len(matches)==1 and matches[0].isfile()
  stream=tar.extractfile(matches[0]);assert stream is not None
  with zstandard.ZstdDecompressor().stream_reader(stream) as reader,payload.open('wb') as out:shutil.copyfileobj(reader,out,8*1024*1024)
 result={'passed':True,'member':'pred.h5ad.zst','full_CSR_comparison':compare_h5(staged,payload),'scope':'Every prediction value/index/pointer and canonical official cell indices and all gene/context/target axis labels in the compressed package match the losslessly staged source. Official validate_vcc checks wrapper presence only.'}
 Path(receipt_path).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 payload.unlink();return result


