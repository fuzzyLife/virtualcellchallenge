"""Metadata-only bounded range audit of author CD4 donor/condition power."""
from pathlib import Path
import io,json,hashlib,requests
from concurrent.futures import ThreadPoolExecutor
import h5py,numpy as np,pandas as pd
pd.options.future.infer_string=False
from prepare_reference_controls import column
import build_submission as b
ROOT=b.ROOT;A=ROOT/'experiments/v22_state_response/analysis';R=ROOT/'reference/v22_cd4'
class Remote(io.RawIOBase):
 def __init__(self,item):
  self.item=item;self.url=item['links'].get('self') or item['links']['content'];self.size=item['size'];self.position=0;self.blocks={};self.block_size=1024**2;self.used=0;self.folder=R/(item['key']+'_metadata');self.folder.mkdir(exist_ok=True)
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.position
 def seek(self,n,w=0):self.position=n if w==0 else self.position+n if w==1 else self.size+n;return self.position
 def read(self,n=-1):
  end=min(self.size,self.position+n) if n>=0 else self.size;pieces=[]
  while self.position<end:
   bi=self.position//self.block_size;lo=bi*self.block_size;hi=min(self.size,lo+self.block_size)
   if bi not in self.blocks:
    p=self.folder/f'{bi}.bin'
    if p.exists():data=p.read_bytes()
    else:
     r=requests.get(self.url,headers={'Range':f'bytes={lo}-{hi-1}'},timeout=(20,60));r.raise_for_status();assert r.status_code==206 and r.headers['Content-Range']==f'bytes {lo}-{hi-1}/{self.size}';data=r.content;p.write_bytes(data)
    assert len(data)==hi-lo;self.used+=len(data)
    if self.used>48*1024**2:raise RuntimeError('Metadata read budget exceeded')
    self.blocks[bi]=data
   stop=min(end,hi);pieces.append(self.blocks[bi][self.position-lo:stop-lo]);self.position=stop
  return b''.join(pieces)
 def readinto(self,buf):v=self.read(len(buf));buf[:len(v)]=v;return len(v)
def inspect(item):
 rr=Remote(item)
 with h5py.File(rr,'r') as h:
  values=column(h['obs']['perturbation']);targets=pd.read_csv(ROOT/'pert_counts.csv').target_gene.to_numpy(str);counts=pd.Series(values).value_counts();selected={t:int(counts.get(t,0)) for t in targets};row={'file':item['key'],'source_bytes':item['size'],'expected_md5':item['checksum'],'overlap':sum(v>0 for v in selected.values()),'targets400':sum(v>=400 for v in selected.values()),'targets200':sum(v>=200 for v in selected.values()),'counts':selected,'metadata_range_bytes':rr.used,'source_blocks_sha256':{str(k):hashlib.sha256(v).hexdigest() for k,v in rr.blocks.items()},'treatment_expression_read':False}
 (A/(item['key']+'_power.json')).write_text(json.dumps(row,indent=2),encoding='utf-8');print(item['key'],'targets400',row['targets400'],'targets200',row['targets200'],flush=True);return row
if __name__=='__main__':
 record=json.loads((A/'cipher_zenodo_record.json').read_text());names=['Marson2025_D2_Rest_filtered.h5ad','Marson2025_D1_Stim48hr_filtered.h5ad','Marson2025_D4_Rest_filtered.h5ad'];items=[v for v in record['files'] if v['key'] in names]
 with ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(inspect,items))
 (A/'additional_cd4_power.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
