"""Fetch metadata-selected fresh donor counts in bounded resumable public ranges."""
from pathlib import Path
import json,hashlib,requests,shutil,time
from concurrent.futures import ThreadPoolExecutor
import build_submission as b
ROOT=b.ROOT;A=ROOT/'experiments/v25_cross_domain/analysis';R=ROOT/'reference/v25_sources'
def main():
 record=json.loads((ROOT/'experiments/v24_cd4_response/analysis/cipher_zenodo_record.json').read_text());item=next(v for v in record['files'] if v['key']=='Marson2025_D3_Rest_filtered.h5ad');size=item['size'];url=item['links']['self'];dest=R/item['key'];folder=R/'D3_Rest_parts';folder.mkdir(exist_ok=True);choice=A/'reserve_source_choice.json'
 if not choice.exists():b.write_json(choice,{'source':item,'metadata_power_sha256':b.sha256(A/(item['key']+'_power.json')),'reason':'14 target groups>=400 by metadata; new D3 raw donor. Reserve excludes whole CD4 family, including all related public summary and raw sources. No treatment expression inspected in source choice.'})
 def fetch(task):
  lo,hi=task;path=folder/f'{lo}_{hi}.bin'
  for attempt in range(8):
   n=path.stat().st_size if path.exists() else 0;assert n<=hi-lo
   if n==hi-lo:return path
   try:
    with requests.get(url,headers={'Range':f'bytes={lo+n}-{hi-1}'},stream=True,timeout=(15,30)) as r:
     r.raise_for_status();assert r.status_code==206 and r.headers['Content-Range']==f'bytes {lo+n}-{hi-1}/{size}'
     with path.open('ab' if n else 'wb') as f:
      for block in r.iter_content(256*1024):
       if block:f.write(block)
    assert path.stat().st_size==hi-lo;return path
   except (requests.RequestException,AssertionError):
    print('V25 source retry',lo,attempt+1,flush=True)
    if attempt==7:raise
    time.sleep(min(attempt+1,5))
 if not dest.exists():
  tasks=[(lo,min(lo+16*1024**2,size)) for lo in range(0,size,16*1024**2)]
  with ThreadPoolExecutor(max_workers=16) as pool:
   for i,_ in enumerate(pool.map(fetch,tasks)):
    if (i+1)%16==0:print('V25 new donor ranges',i+1,'/',len(tasks),flush=True)
  pending=dest.with_suffix('.h5ad.partial')
  with pending.open('wb') as f:
   for lo,hi in tasks:
    with (folder/f'{lo}_{hi}.bin').open('rb') as src:shutil.copyfileobj(src,f,8*1024**2)
  assert pending.stat().st_size==size
  with pending.open('rb') as f:assert hashlib.file_digest(f,'md5').hexdigest()==item['checksum'].split(':')[1]
  pending.replace(dest)
 else:
  with dest.open('rb') as f:assert hashlib.file_digest(f,'md5').hexdigest()==item['checksum'].split(':')[1]
 b.write_json(A/'holdout_download.json',{'file':item,'bytes':size,'md5_verified':True,'source_choice_sha256':b.sha256(choice),'collector_sha256':b.sha256(Path(__file__))});print('V25 new donor full MD5 verified',flush=True)
if __name__=='__main__':main()
