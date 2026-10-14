"""Resume the same MD5-pinned raw holdout using smaller bounded parallel ranges."""
from pathlib import Path
import json,hashlib,requests,shutil,time
from concurrent.futures import ThreadPoolExecutor
import build_submission as b
ROOT=b.ROOT;A=ROOT/'experiments/v24_cd4_response/analysis';R=ROOT/'reference/v24_cd4'
def main():
 rec=json.loads((A/'cipher_zenodo_record.json').read_text());item=next(v for v in rec['files'] if v['key']=='Marson2025_D2_Stim48hr_filtered.h5ad');size=item['size'];url=item['links']['self'];dest=R/item['key'];chunks=R/'bounded_download_chunks';chunks.mkdir(exist_ok=True);amend=A/'download_transport_amendment.json'
 if not amend.exists():
  prefix=[(R/'D2_Stim48_download_parts'/f'{i}.bin').stat().st_size for i in range(4)];b.write_json(amend,{'preserved_prefix_bytes':prefix,'parts':4,'workers':16,'range_bytes':16*1024**2,'reason':'Large individual streaming requests stalled; transport-only resume preserves each exact existing prefix and author full-file MD5. No outcome-dependent data/model change.'})
 prefix=json.loads(amend.read_text())['preserved_prefix_bytes'];tasks=[]
 for i in range(4):
  lo=size*i//4;hi=size*(i+1)//4;path=R/'D2_Stim48_download_parts'/f'{i}.bin';assert path.stat().st_size==prefix[i]
  for start in range(lo+prefix[i],hi,16*1024**2):tasks.append((start,min(start+16*1024**2,hi)))
 def fetch(task):
  lo,hi=task;file=chunks/f'{lo}_{hi}.bin'
  for attempt in range(8):
   n=file.stat().st_size if file.exists() else 0
   assert n<=hi-lo
   if n==hi-lo:return file
   try:
    with requests.get(url,headers={'Range':f'bytes={lo+n}-{hi-1}'},stream=True,timeout=(15,30)) as response:
     response.raise_for_status();assert response.status_code==206 and response.headers['Content-Range']==f'bytes {lo+n}-{hi-1}/{size}'
     with file.open('ab' if n else 'wb') as out:
      for block in response.iter_content(256*1024):
       if block:out.write(block)
    assert file.stat().st_size==hi-lo;print('V24 bounded range complete',lo,hi,flush=True);return file
   except (requests.RequestException,AssertionError) as error:
    print('V24 range retry',lo,attempt+1,type(error).__name__,flush=True)
    if attempt==7:raise
    time.sleep(min(attempt+1,5))
 if not dest.exists():
  with ThreadPoolExecutor(max_workers=16) as pool:list(pool.map(fetch,tasks))
  pending=dest.with_suffix('.h5ad.partial')
  with pending.open('wb') as out:
   for i in range(4):
    with (R/'D2_Stim48_download_parts'/f'{i}.bin').open('rb') as f:shutil.copyfileobj(f,out,8*1024**2)
    lo=size*i//4;hi=size*(i+1)//4
    for start in range(lo+prefix[i],hi,16*1024**2):
     end=min(start+16*1024**2,hi)
     with (chunks/f'{start}_{end}.bin').open('rb') as f:shutil.copyfileobj(f,out,8*1024**2)
  assert pending.stat().st_size==size
  with pending.open('rb') as f:assert hashlib.file_digest(f,'md5').hexdigest()==item['checksum'].split(':')[1]
  pending.replace(dest)
 with dest.open('rb') as f:assert hashlib.file_digest(f,'md5').hexdigest()==item['checksum'].split(':')[1]
 b.write_json(A/'holdout_download.json',{'file':item,'bytes':size,'md5_verified':True,'metadata_only_selection':'Untouched D2 Stim48hr has13 competition targets>=400; all condition response rows excluded from validation training.','download_script_sha256':b.sha256(Path(__file__)),'transport_amendment_sha256':b.sha256(amend)});print('V24 full holdout MD5 verified',flush=True)
if __name__=='__main__':main()
