"""Extract competition-target DE rows via ETag-pinned HTTP byte ranges."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib, io, json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
ROOT=Path(__file__).resolve().parent
BASE='https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/'
OBJECTS={'merged':('GWCD4i.DE_stats.h5ad',16786240107,'"c9ff52fcc6d6ce8a387a76dc757a5b97-2002"'),
 'donors':('GWCD4i.DE_stats.by_donors.h5mu',16866278447,'"1f405d8c2d202c6e845c1d11e59f2e47-2011"')}
def session():
 s=requests.Session();s.mount('https://',HTTPAdapter(max_retries=Retry(total=3,backoff_factor=.5,status_forcelist=[429,500,502,503,504])));return s
def get_range(client,source,start,stop):
 name,size,etag=source
 if not 0<=start<stop<=size:raise ValueError('Invalid byte range')
 with client.get(BASE+name,headers={'Range':f'bytes={start}-{stop-1}','If-Match':etag},timeout=(15,90),stream=True) as r:
  if r.status_code!=206 or r.headers.get('Content-Range')!=f'bytes {start}-{stop-1}/{size}':raise RuntimeError(f'Range rejected: {r.status_code}')
  data=r.content
 if len(data)!=stop-start:raise RuntimeError('Incomplete HTTP range')
 return data
class MetadataReader(io.RawIOBase):
 """Seekable remote HDF5 metadata reader, limited to 64 MiB."""
 def __init__(self,source):self.source=source;self.position=0;self.cache={};self.downloaded=0;self.client=session()
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.position
 def seek(self,offset,whence=0):
  self.position=offset if whence==0 else self.position+offset if whence==1 else self.source[1]+offset;return self.position
 def read(self,count=-1):
  stop=min(self.source[1],self.position+count) if count>=0 else self.source[1];parts=[]
  while self.position<stop:
   block=self.position//262144;lo=block*262144;hi=min(self.source[1],(block+1)*262144)
   if block not in self.cache:
    if self.downloaded+hi-lo>64*1024**2:raise RuntimeError('Metadata budget exceeded')
    self.cache[block]=get_range(self.client,self.source,lo,hi);self.downloaded+=hi-lo
   end=min(stop,hi);parts.append(self.cache[block][self.position-lo:end-lo]);self.position=end
  return b''.join(parts)
 def readinto(self,b):
  data=self.read(len(b));b[:len(data)]=data;return len(data)
def column(obj):
 if isinstance(obj,h5py.Group):
  values=column(obj['categories']);codes=obj['codes'][:];return np.where(codes>=0,values[np.maximum(codes,0)],'')
 return obj.asstr()[:] if obj.dtype.kind in 'OSU' else obj[:]
def extract(group,source,targets,destination):
 metadata=pd.DataFrame({key:column(group['obs'][key]) for key in group['obs']})
 selected=np.flatnonzero(metadata.target_contrast_gene_name.isin(targets));metadata=metadata.iloc[selected].copy();metadata['source_row']=selected
 symbols=column(group['var/gene_name']);runs=[]
 for row in selected:
  if runs and row-runs[-1][1]<=3:runs[-1][1]=int(row)+1
  else:runs.append([int(row),int(row)+1])
 arrays={'genes':np.asarray(symbols,dtype=str),'source_rows':selected};planned=0;downloaded=0
 for layer in ['log_fc','adj_p_value','lfcSE']:
  ds=group['layers'][layer];offset=ds.id.get_offset();dtype=ds.dtype;shape=ds.shape
  if offset is None or ds.chunks is not None or dtype.kind!='f':raise RuntimeError('Expected contiguous floating DE layer')
  width=shape[1]*dtype.itemsize;planned+=sum((end-start)*width for start,end in runs)
  if planned>400*1024**2:raise RuntimeError('Per-modality budget exceeded')
  def fetch(run):
   start,end=run
   with session() as client:raw=get_range(client,source,offset+start*width,offset+end*width)
   values=np.frombuffer(raw,dtype=dtype).reshape(end-start,shape[1]);positions=selected[(selected>=start)&(selected<end)]-start
   return values[positions].copy(),len(raw)
  with ThreadPoolExecutor(max_workers=6) as pool:results=list(pool.map(fetch,runs))
  arrays[layer]=np.concatenate([value for value,_ in results]).astype(np.float32);downloaded+=sum(count for _,count in results)
  print(f'{destination.stem}: {layer}, {len(selected)} rows, {downloaded/1e6:.1f} MB',flush=True)
 np.savez_compressed(destination,**arrays);metadata.to_csv(destination.with_suffix('.csv'),index=False)
 return {'subset':destination.name,'rows':len(selected),'genes':len(symbols),'bytes_downloaded':downloaded,'sha256':hashlib.sha256(destination.read_bytes()).hexdigest(),'source':BASE+source[0],'source_size':source[1],'source_etag':source[2]}
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--inspect',action='store_true');args=parser.parse_args()
 directory=ROOT/'reference/cd4_probe';directory.mkdir(parents=True,exist_ok=True);targets=pd.read_csv(ROOT/'pert_counts.csv').target_gene.astype(str).to_numpy();report=[]
 for kind,source in OBJECTS.items():
  reader=MetadataReader(source)
  with h5py.File(reader,'r') as f:
   groups={'merged':f} if kind=='merged' else dict(f['mod'].items());print(kind,list(groups),flush=True)
   if args.inspect:
    for name,group in groups.items():print(name,{key:(ds.shape,str(ds.dtype),ds.id.get_offset()) for key,ds in group['layers'].items()},flush=True)
    continue
   if kind=='donors':
    first=sorted(groups)[0];members=set(first.split('_'));other=[name for name in groups if not members.intersection(name.split('_'))]
    if len(other)!=1:raise RuntimeError('Complementary donor pairs not found')
    groups={name:groups[name] for name in [first,other[0]]}
   for name,group in groups.items():report.append(extract(group,source,targets,directory/f'{name}.npz'))
  reader.close()
 if not args.inspect:(ROOT/'analysis/model_audit/cd4_extraction.json').write_text(json.dumps(report,indent=2))
if __name__=='__main__':main()