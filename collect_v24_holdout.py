"""Download the metadata-selected powered CD4 Stim48 reserve, with range/MD5 checks."""
from pathlib import Path
import json,hashlib,requests,shutil,time
from concurrent.futures import ThreadPoolExecutor
import build_submission as b
ROOT=b.ROOT;A=ROOT/'experiments/v24_cd4_response/analysis';R=ROOT/'reference/v24_cd4'
def main():
 R.mkdir(parents=True,exist_ok=True);record=json.loads((A/'cipher_zenodo_record.json').read_text());item=next(v for v in record['files'] if v['key']=='Marson2025_D2_Stim48hr_filtered.h5ad');dest=R/item['key'];folder=R/'D2_Stim48_download_parts';folder.mkdir(exist_ok=True);url=item['links'].get('self') or item['links']['content'];size=item['size'];parts=4
 def fetch(i):
  lo=size*i//parts;hi=size*(i+1)//parts;path=folder/f'{i}.bin';start=path.stat().st_size if path.exists() else 0;assert start<=hi-lo
  while start<hi-lo:
   end=min(lo+start+32*1024**2,hi)-1
   for attempt in range(6):
    try:
     with requests.get(url,headers={'Range':f'bytes={lo+start}-{end}'},stream=True,timeout=(15,30)) as response:
      response.raise_for_status();assert response.status_code==206 and response.headers['Content-Range']==f'bytes {lo+start}-{end}/{size}'
      with path.open('ab' if start else 'wb') as f:
       for block in response.iter_content(1024**2):
        if block:f.write(block)
     start=path.stat().st_size;assert start<=hi-lo
     print('CD4 Stim48 part',i,start,'/',hi-lo,flush=True);break
    except (requests.RequestException,AssertionError) as error:
     start=path.stat().st_size if path.exists() else 0
     print('CD4 download bounded retry',i,attempt+1,type(error).__name__,str(error)[:180],flush=True)
     if attempt==5:raise
     time.sleep(min(2*(attempt+1),10));end=min(lo+start+32*1024**2,hi)-1
  assert path.stat().st_size==hi-lo;return path
 if not dest.exists():
  with ThreadPoolExecutor(max_workers=parts) as pool:paths=list(pool.map(fetch,range(parts)))
  pending=dest.with_suffix('.h5ad.partial')
  with pending.open('wb') as out:
   for path in paths:
    with path.open('rb') as src:shutil.copyfileobj(src,out,8*1024**2)
  assert pending.stat().st_size==size
  with pending.open('rb') as f:assert hashlib.file_digest(f,'md5').hexdigest()==item['checksum'].split(':')[1]
  pending.replace(dest)
 with dest.open('rb') as f:assert hashlib.file_digest(f,'md5').hexdigest()==item['checksum'].split(':')[1]
 (A/'holdout_download.json').write_text(json.dumps({'file':item,'bytes':size,'md5_verified':True,'metadata_only_selection':'Untouched D2 Stim48hr has13 competition targets>=400 by metadata; its treatment response rows are excluded from the new expert and baseline fitting.','download_script_sha256':b.sha256(Path(__file__))},indent=2),encoding='utf-8');print('CD4 Stim48 raw source MD5 verified',flush=True)
if __name__=='__main__':main()
