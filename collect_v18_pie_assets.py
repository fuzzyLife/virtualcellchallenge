"""Fetch pinned PIE checkpoint/dense sources and only the 300 needed NCBI token slices."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,hashlib,io,json,time
import numpy as np,requests
ROOT=Path(__file__).resolve().parent;R=ROOT/'reference/v18_pie';A=ROOT/'experiments/v18_pie_expert/analysis'
REV='cb1aaa4e7655605bdc70a9bd77bbd62016b8c7d7'
MODEL_REV='055c7a2cabdaff5659121c64fa2b609c231fc5be'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(8*1024**2),b''):h.update(x)
 return h.hexdigest()
def write(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')
def get_range(url,lo,hi):
 for attempt in range(4):
  try:
   with requests.get(url+('?slice='+str(lo) if '?' not in url else '&slice='+str(lo)),headers={'Range':f'bytes={lo}-{hi}'},stream=True,timeout=(30,180)) as r:
    r.raise_for_status();assert r.status_code==206,(r.status_code,dict(r.headers));assert r.headers['Content-Range'].startswith(f'bytes {lo}-{hi}/'),r.headers.get('Content-Range')
    data=bytearray()
    for v in r.iter_content(1024**2):
     data.extend(v)
     if len(data)>hi-lo+1:raise ValueError('Server ignored bounded range')
    assert len(data)==hi-lo+1;return bytes(data)
  except Exception:
   if attempt==3:raise
   time.sleep(2**attempt)
def fetch(url,path,size,expected=None):
 path.parent.mkdir(parents=True,exist_ok=True)
 if path.exists():
  assert path.stat().st_size==size and (expected is None or sha(path)==expected);return
 tmp=path.with_suffix(path.suffix+'.partial')
 with requests.get(url,stream=True,timeout=(30,180)) as r:
  r.raise_for_status()
  with tmp.open('wb') as f:
   n=0
   for chunk in r.iter_content(4*1024**2):
    f.write(chunk);n+=len(chunk)
    if n>size:raise ValueError('Unexpected download size')
 assert tmp.stat().st_size==size and (expected is None or sha(tmp)==expected)
 tmp.replace(path);print('Verified',path.relative_to(ROOT),size,flush=True)
def source_url(path):return f'https://huggingface.co/datasets/arcinstitute/PIE_sources/resolve/{REV}/{path}'
def main():
 A.mkdir(parents=True,exist_ok=True);files=json.loads((ROOT/'analysis/post_v17_research/PIE_sources_files.json').read_text());files={v['path']:v for v in files if v['type']=='file'}
 dense=['esm2','string_space','depmap_gene_effect','context_text','perturbation_text','gene_text'];todo=[]
 for name in dense+['ncbi_text']:
  for file in (['meta.json','offsets.npy'] if name=='ncbi_text' else ['meta.json','embeddings.npy']):
   rel=name+'/'+file;rec=files[rel];todo.append((source_url(rel),R/'sources'/rel,rec['size'],rec.get('lfs',{}).get('oid')))
 todo.append((f'https://huggingface.co/arcinstitute/PIE_replogle_wdataset/resolve/{MODEL_REV}/k562/best_auprc.ckpt',R/'best_auprc.ckpt',949595187,'b9d5a0b3943666d892a56e7a9368f7850826b3f1fe6879e63a1a590f8add7d84'))
 with ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(lambda v:fetch(*v),todo))
 # Byte ranges retain every original token for each selected target; no truncation/zero approximation.
 meta=json.loads((R/'sources/ncbi_text/meta.json').read_text());off=np.load(R/'sources/ncbi_text/offsets.npy',allow_pickle=False)
 import pandas as pd
 targets=pd.read_csv(ROOT/'pert_counts.csv').iloc[:,0].astype(str).tolist();keys=[k for k in meta['keys'] if k in targets];assert len(keys)>0
 url=source_url('ncbi_text/embeddings.npy');header=get_range(url,0,127);buf=io.BytesIO(header);ver=np.lib.format.read_magic(buf);shape,fortran,dtype=(np.lib.format.read_array_header_1_0(buf) if ver==(1,0) else np.lib.format.read_array_header_2_0(buf));start=buf.tell();assert not fortran and shape[1]==meta['dim'] and dtype==np.dtype(meta['dtype']);stride=shape[1]*dtype.itemsize
 chunks=R/'ncbi_ranges';chunks.mkdir(exist_ok=True);index={k:i for i,k in enumerate(meta['keys'])}
 def one(k):
  i=index[k];lo=start+int(off[i])*stride;hi=start+int(off[i+1])*stride-1;p=chunks/(k+'.bin');receipt=p.with_suffix('.json')
  if p.exists():
   rec=json.loads(receipt.read_text());assert rec['sha256']==sha(p) and p.stat().st_size==hi-lo+1
  else:
   data=get_range(url,lo,hi);p.write_bytes(data);write(receipt,{'key':k,'lo':lo,'hi':hi,'url':url,'sha256':sha(p),'tokens':int(off[i+1]-off[i])})
  return json.loads(receipt.read_text())
 with ThreadPoolExecutor(max_workers=6) as pool:receipts=list(pool.map(one,keys))
 compact=R/'sources/ncbi_subset';compact.mkdir(exist_ok=True);offsets=np.r_[0,np.cumsum([v['tokens'] for v in receipts])].astype(np.int64)
 out=np.lib.format.open_memmap(compact/'embeddings.npy',mode='w+',dtype=dtype,shape=(int(offsets[-1]),shape[1]))
 for i,k in enumerate(keys):out[offsets[i]:offsets[i+1]]=np.fromfile(chunks/(k+'.bin'),dtype=dtype).reshape(-1,shape[1])
 out.flush();del out;np.save(compact/'offsets.npy',offsets);meta['keys']=keys;meta['provenance']['v18_subset']={'source_revision':REV,'range_receipts':'experiments/v18_pie_expert/analysis/ncbi_ranges.json','no_tokens_removed':True};write(compact/'meta.json',meta)
 write(A/'ncbi_ranges.json',{'source_array_bytes':files['ncbi_text/embeddings.npy']['size'],'source_lfs_sha256':files['ncbi_text/embeddings.npy']['lfs']['oid'],'npy_header_hex':header.hex(),'ranges':receipts,'total_payload_bytes':sum(v['hi']-v['lo']+1 for v in receipts)})
 records=[]
 for p in sorted((R/'sources').rglob('*')):
  if p.is_file():records.append({'path':str(p.relative_to(ROOT)).replace('\\','/'),'bytes':p.stat().st_size,'sha256':sha(p)})
 for p in sorted((R/'code/src/pie').rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts:records.append({'path':str(p.relative_to(ROOT)).replace('\\','/'),'bytes':p.stat().st_size,'sha256':sha(p)})
 records.append({'path':'reference/v18_pie/best_auprc.ckpt','bytes':949595187,'sha256':sha(R/'best_auprc.ckpt')})
 write(A/'assets.json',{'source_revision':REV,'model_revision':MODEL_REV,'code_revision':'21a526add9ff3a81b4f962d9b16bbae570052c63','files':records});print('PIE assets complete; NCBI bytes',sum(v['hi']-v['lo']+1 for v in receipts),flush=True)
if __name__=='__main__':main()
