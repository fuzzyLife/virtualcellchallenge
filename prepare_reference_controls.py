"""Bounded public control-profile extraction; never substitute DE for expression."""
from pathlib import Path
import argparse, hashlib, io, json
import h5py
import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
ROOT=Path(__file__).resolve().parent
SOURCES={
'H1':'https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/train/adata_Training.h5ad?generation=1765904883947296',
'CD4':'https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.pseudobulk_merged.h5ad'}
def column(obj):
 if isinstance(obj,h5py.Group):
  if 'categories' in obj:
   categories=column(obj['categories']);codes=obj['codes'][:]
   return np.where(codes>=0,categories[np.maximum(codes,0)],'')
  if 'values' in obj:
   values=column(obj['values']).copy();values[obj['mask'][:]]='';return values
  raise ValueError(f'Unsupported column {obj.name}')
 return obj.asstr()[:] if obj.dtype.kind in 'OSU' else obj[:]
def index(group):
 name=group.attrs.get('_index','_index')
 if isinstance(name,bytes):name=name.decode()
 return column(group[name])
class RangeReader(io.RawIOBase):
 """Persist SHA-verified blocks, pin source ETag, and bound bytes used."""
 def __init__(self,url,cache,budget=1024**3):
  self.url=url;self.cache=Path(cache);self.cache.mkdir(parents=True,exist_ok=True)
  self.client=requests.Session();self.client.mount('https://',HTTPAdapter(max_retries=Retry(total=3,backoff_factor=.5,status_forcelist=[429,500,502,503,504])))
  with self.client.get(url,headers={'Range':'bytes=0-0'},timeout=(15,90),stream=True) as r:
   if r.status_code!=206:raise RuntimeError(f'Range unavailable: {r.status_code}')
   self.size=int(r.headers['Content-Range'].split('/')[-1]);self.etag=r.headers['ETag']
  self.identity={'url':url,'bytes':self.size,'etag':self.etag};record=self.cache/'identity.json'
  if record.exists() and json.loads(record.read_text())!=self.identity:raise RuntimeError('Remote source changed')
  record.write_text(json.dumps(self.identity,indent=2),encoding='utf-8')
  self.position=0;self.block_size=262144;self.budget=budget;self.blocks={};self.used=0;self.downloaded=0
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.position
 def seek(self,offset,whence=0):
  self.position=offset if whence==0 else self.position+offset if whence==1 else self.size+offset;return self.position
 def read(self,count=-1):
  end=min(self.size,self.position+count) if count>=0 else self.size;parts=[]
  while self.position<end:
   block=self.position//self.block_size;lo=block*self.block_size;hi=min(self.size,lo+self.block_size)
   if block not in self.blocks:
    if self.used+hi-lo>self.budget:raise RuntimeError('Remote byte budget exceeded')
    path=self.cache/f'{block:09d}.bin';checksum=path.with_suffix('.sha256')
    if path.exists() and checksum.exists():
     raw=path.read_bytes()
     if len(raw)!=hi-lo or hashlib.sha256(raw).hexdigest()!=checksum.read_text().strip():raise RuntimeError('Cached range checksum mismatch')
    else:
     with self.client.get(self.url,headers={'Range':f'bytes={lo}-{hi-1}','If-Match':self.etag},timeout=(15,90)) as r:
      if r.status_code!=206 or r.headers.get('Content-Range')!=f'bytes {lo}-{hi-1}/{self.size}':raise RuntimeError('Range response mismatch')
      raw=r.content
     if len(raw)!=hi-lo:raise RuntimeError('Truncated range')
     path.write_bytes(raw);checksum.write_text(hashlib.sha256(raw).hexdigest());self.downloaded+=len(raw)
    self.used+=len(raw);self.blocks[block]=raw
    if len(self.blocks)%128==0:print(f'Ranges used {self.used/1e6:.1f} MB',flush=True)
   stop=min(end,hi);parts.append(self.blocks[block][self.position-lo:stop-lo]);self.position=stop
  return b''.join(parts)
 def readinto(self,b):
  value=self.read(len(b));b[:len(value)]=value;return len(value)
 def close(self):self.client.close();super().close()
def inspect(source):
 reader=RangeReader(SOURCES[source],ROOT/f'reference/control_profiles/ranges_{source}',budget=128*1024**2)
 with h5py.File(reader,'r') as f:
  print(source,'keys',list(f),flush=True)
  for axis in ['obs','var']:
   print(axis,dict(f[axis].attrs),flush=True)
   for name,obj in f[axis].items():
    values=column(obj);print(name,values.shape,values[:3].tolist(),'sample values',pd.Series(values).value_counts().head(6).to_dict(),flush=True)
  def show(name,obj):
   if isinstance(obj,h5py.Dataset):print(name,obj.shape,str(obj.dtype),obj.chunks,obj.compression,flush=True)
  for name in ['X','layers']:
   if name in f:
    print(name,dict(f[name].attrs),flush=True)
    if isinstance(f[name],h5py.Group):
     for key in f[name]:show(key,f[name][key])
    else:show(name,f[name])
 print('remote bytes',reader.used,flush=True);reader.close()
def save_profile(name,genes,counts,record,output):
 if len(genes)!=len(counts) or not np.isfinite(counts).all() or (counts<0).any() or counts.sum()<=0:raise ValueError('Invalid reference expression')
 frame=pd.DataFrame({'gene_name':np.asarray(genes,dtype=str),'control_count_sum':counts})
 frame=frame.groupby('gene_name',sort=False,as_index=False).control_count_sum.sum()
 frame['pooled_cpm']=1e6*frame.control_count_sum/frame.control_count_sum.sum()
 path=output/f'{name}.csv';frame.to_csv(path,index=False)
 record.update({'name':name,'genes':len(frame),'profile_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'normalization':'pool original counts, normalize to CPM; sum duplicate gene symbols'})
 (output/f'{name}.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
 print(name,record,flush=True)

def prepare_local(output):
 import pyarrow.parquet as pq
 import pilot_orion
 # Verified original K562 pseudobulk means, weighted by source cell membership.
 path=ROOT/'reference/K562_gwps_raw_bulk_01.h5ad'
 if hashlib.md5(path.read_bytes()).hexdigest()!='4570b53c9d62ff6df281e622f0350060':raise RuntimeError('Pinned K562 checksum mismatch')
 with h5py.File(path,'r') as f:
  labels=f['obs/gene_transcript'].asstr()[:];n=f['obs/num_cells_filtered'][:]
  targets=np.array([s.split('_',2)[1] if len(s.split('_',2))==3 and s.split('_',1)[0].isdigit() else '' for s in labels])
  rows=np.flatnonzero((targets=='non-targeting')&f['obs/core_control'][:].astype(bool)&(n>=30))
  ds=f['var/gene_name'];genes=f[ds.attrs['categories']].asstr()[:][ds[:]]
  if not len(rows):raise ValueError('No K562 controls')
  counts=np.sum(f['X'][rows,:]*n[rows,None],axis=0)
  save_profile('K562',genes,counts,{'source':json.loads((ROOT/'reference/provenance.json').read_text()),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'control_groups':len(rows),'control_cells':int(n[rows].sum()),'selection':'non-targeting, core_control, at least 30 cells; raw pseudobulk mean times membership'},output)
  pd.DataFrame({'source_row':rows,'label':labels[rows],'n_cells':n[rows]}).to_csv(output/'K562_selected_controls.csv',index=False)
 metadata_path=ROOT/'reference/orion_pilot/gene_metadata.parquet'
 if hashlib.sha256(metadata_path.read_bytes()).hexdigest()!=pilot_orion.FILES['metadata/gene_metadata.parquet'][1]:raise RuntimeError('Pinned Orion gene metadata checksum mismatch')
 meta=pq.read_table(metadata_path).to_pandas()
 for name in ['HCT116','HEK293T']:
  path=ROOT/f'reference/orion_pilot/{name}_Batch1.parquet'
  expected=pilot_orion.FILES['data/'+path.name][1]
  if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:raise RuntimeError('Pinned Orion checksum mismatch')
  table=pq.read_table(path,columns=['gene_target','pass_guide_filter','total_counts','gene_token_id','gene_expression'])
  labels=np.array(table['gene_target'].to_pylist());n=table['total_counts'].to_numpy()
  rows=np.flatnonzero(np.array([str(x).lower().startswith('non-targeting') or str(x).lower() in ['ntc','control'] for x in labels])&(table['pass_guide_filter'].to_numpy()>0)&(n>0))
  counts=np.zeros(int(meta.gene_token_id.max())+1)
  for row in rows:
   tokens=np.array(table['gene_token_id'][int(row)].as_py(),dtype=int);values=np.array(table['gene_expression'][int(row)].as_py(),dtype=float)
   if not np.isfinite(values).all() or (values<0).any() or (values!=np.floor(values)).any() or values.sum()!=n[row]:raise ValueError('Invalid Orion counts')
   np.add.at(counts,tokens,values)
  save_profile(name,meta.gene_name.to_numpy(),counts[meta.gene_token_id.to_numpy()],{'url':f'https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion/resolve/{pilot_orion.REVISION}/data/{path.name}','revision':pilot_orion.REVISION,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'gene_metadata_sha256':hashlib.sha256((ROOT/'reference/orion_pilot/gene_metadata.parquet').read_bytes()).hexdigest(),'control_cells':len(rows),'scope':'All passing controls in Batch1 only, not the full atlas'},output)
  pd.DataFrame({'source_row':rows,'target_gene':labels[rows],'total_counts':n[rows]}).to_csv(output/f'{name}_selected_controls.csv',index=False)

def select_stratified(meta,mask,keys,limit):
 eligible=meta.loc[mask].copy();eligible['source_row']=np.flatnonzero(mask)
 return eligible.groupby(keys,sort=True,observed=True,group_keys=False).head(limit).sort_values('source_row')

def prefetch_sparse(reader,x,indptr,rows):
 """Plan only the HDF5 storage chunks covering selected control rows."""
 from concurrent.futures import ThreadPoolExecutor
 blocks=set(reader.blocks)
 for key in ['data','indices']:
  ds=x[key];print('Storage',key,ds.shape,ds.chunks,ds.compression,flush=True)
  spans=[]
  if ds.chunks is None:
   offset=ds.id.get_offset()
   if offset is None:raise ValueError('Unknown sparse storage')
   spans=[(offset+int(indptr[r])*ds.dtype.itemsize,offset+int(indptr[r+1])*ds.dtype.itemsize) for r in rows]
  else:
   chunks=set()
   for r in rows:
    lo,hi=int(indptr[r]),int(indptr[r+1])
    chunks.update(range(lo//ds.chunks[0],(hi-1)//ds.chunks[0]+1))
   for chunk in sorted(chunks):
    info=ds.id.get_chunk_info_by_coord((chunk*ds.chunks[0],))
    if info.byte_offset is None:raise ValueError('Unallocated selected sparse chunk')
    spans.append((info.byte_offset,info.byte_offset+info.size))
  for lo,hi in spans:
   blocks.update(range(lo//reader.block_size,(hi-1)//reader.block_size+1))
 planned=sum(min(reader.block_size,reader.size-b*reader.block_size) for b in blocks)
 if planned>reader.budget:raise RuntimeError(f'Planned source ranges {planned} exceed byte budget')
 missing=[b for b in sorted(blocks) if not (reader.cache/f'{b:09d}.bin').exists()]
 print('Prefetch',len(missing),'blocks; bounded planned bytes',planned,flush=True)
 def fetch(block):
  lo=block*reader.block_size;hi=min(reader.size,lo+reader.block_size)
  # A separate retrying session per worker request avoids shared session mutation.
  client=requests.Session();client.mount('https://',HTTPAdapter(max_retries=Retry(total=3,backoff_factor=.5,status_forcelist=[429,500,502,503,504])))
  with client.get(reader.url,headers={'Range':f'bytes={lo}-{hi-1}','If-Match':reader.etag},timeout=(15,90)) as response:
   if response.status_code!=206 or response.headers.get('Content-Range')!=f'bytes {lo}-{hi-1}/{reader.size}':raise RuntimeError('Prefetch range mismatch')
   raw=response.content
  client.close()
  if len(raw)!=hi-lo:raise RuntimeError('Truncated prefetch')
  path=reader.cache/f'{block:09d}.bin';path.write_bytes(raw);path.with_suffix('.sha256').write_text(hashlib.sha256(raw).hexdigest())
  return len(raw)
 with ThreadPoolExecutor(max_workers=8) as pool:
  for i,size in enumerate(pool.map(fetch,missing)):
   reader.downloaded+=size
   if (i+1)%128==0:print('Prefetched',i+1,'/',len(missing),'blocks',flush=True)
 reader.prefetch_planned_bytes=planned

def prepare_remote(name,output,per_stratum=32,budget=2*1024**3):
 reader=RangeReader(SOURCES[name],output/f'ranges_{name}',budget=budget)
 expected_size={'H1':15482497461,'CD4':44566657140}[name]
 if reader.size!=expected_size:raise RuntimeError('Pinned source size changed')
 if name=='CD4' and reader.etag!='"010c14e0af0dccbc2524529d28ca517e-5313"':raise RuntimeError('Pinned CD4 source ETag changed')
 with h5py.File(reader,'r') as f:
  if name=='H1':
   meta=pd.DataFrame({k:column(f['obs'][k]) for k in ['target_gene','guide_id','batch']})
   mask=meta.target_gene.to_numpy()=='non-targeting';keys=['batch']
   genes=index(f['var'])
  else:
   names=['guide_type','perturbed_gene_name','n_cells','keep_min_cells','keep_total_counts','culture_condition','donor_id','10xrun_id','guide_id','total_counts']
   meta=pd.DataFrame({k:column(f['obs'][k]) for k in names})
   mask=(meta.perturbed_gene_name.to_numpy()=='NTC')&(meta.guide_type.to_numpy()=='non-targeting')&(meta.n_cells.to_numpy()>0)&meta.keep_min_cells.to_numpy().astype(bool)&meta.keep_total_counts.to_numpy().astype(bool)
   keys=['culture_condition','donor_id','10xrun_id'];genes=column(f['var/gene_name'])
  chosen=select_stratified(meta,mask,keys,per_stratum);rows=chosen.source_row.to_numpy()
  if not len(rows):raise ValueError('No eligible controls')
  print(name,'selected',len(rows),'of',int(mask.sum()),'strata',chosen.groupby(keys,observed=True).ngroups,flush=True)
  chosen.to_csv(output/f'{name}_selected_controls.csv',index=False)
  x=f['X'];indptr=x['indptr'][:];total=np.zeros(len(genes));condition_totals={}
  prefetch_sparse(reader,x,indptr,rows)
  for i,row in enumerate(rows):
   lo,hi=int(indptr[row]),int(indptr[row+1]);values=x['data'][lo:hi].astype(float);columns=x['indices'][lo:hi]
   if not np.isfinite(values).all() or (values<0).any() or (values!=np.floor(values)).any() or values.sum()<=0:raise ValueError('Invalid remote raw controls')
   counts=np.bincount(columns,weights=values,minlength=len(genes));total+=counts
   if name=='CD4':
    if not np.isclose(counts.sum(),float(meta.total_counts.iloc[row]),rtol=0,atol=1):raise ValueError('CD4 count total mismatch')
    condition=str(meta.culture_condition.iloc[row]);condition_totals.setdefault(condition,np.zeros(len(genes)));condition_totals[condition]+=counts
   if (i+1)%32==0:print(name,'controls extracted',i+1,'bytes used',reader.used,flush=True)
  record={'remote':reader.identity,'remote_range_bytes_used':reader.used,'downloaded_this_run':reader.downloaded,'eligible_controls':int(mask.sum()),'selected_control_rows':len(rows),'per_stratum_limit':per_stratum,'stratification':keys,'selection':'first source-ordered eligible controls per stratum; all competition cells are analyzed','control_cells':len(rows) if name=='H1' else int(chosen.n_cells.sum()),'scope':'H1 2025 train release only' if name=='H1' else 'CD4 non-targeting pseudobulks spanning donors, runs and all three culture conditions'}
  save_profile(name,genes,total,record,output)
  for condition,counts in condition_totals.items():
   part=chosen[chosen.culture_condition.astype(str)==condition]
   save_profile('CD4_'+condition,genes,counts,dict(record,scope='CD4 condition '+condition,selected_control_rows=len(part),control_cells=int(part.n_cells.sum())),output)
 reader.close()

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--inspect',choices=['H1','CD4']);parser.add_argument('--source',choices=['all','local','H1','CD4'],default='all')
 parser.add_argument('--controls-per-stratum',type=int,default=32);parser.add_argument('--budget-mib',type=int,default=2048)
 args=parser.parse_args()
 if args.inspect:inspect(args.inspect);return
 if args.controls_per_stratum<1 or args.budget_mib<128:parser.error('Invalid extraction limits')
 output=ROOT/'reference/control_profiles';output.mkdir(parents=True,exist_ok=True)
 if args.source in ['all','local']:prepare_local(output)
 for name in ['H1','CD4']:
  if args.source in ['all',name]:prepare_remote(name,output,args.controls_per_stratum,args.budget_mib*1024**2)
if __name__=='__main__':main()
