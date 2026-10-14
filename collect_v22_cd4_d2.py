"""Download the metadata-selected powered CD4 D2 reserve, with range/MD5 checks."""
from pathlib import Path
import json,hashlib,requests,shutil
from concurrent.futures import ThreadPoolExecutor
import build_submission as b
ROOT=b.ROOT;A=ROOT/'experiments/v22_state_response/analysis';R=ROOT/'reference/v22_cd4'
def main():
 record=json.loads((A/'cipher_zenodo_record.json').read_text());item=next(v for v in record['files'] if v['key']=='Marson2025_D2_Rest_filtered.h5ad');dest=R/item['key'];folder=R/'D2_download_parts';folder.mkdir(exist_ok=True);url=item['links'].get('self') or item['links']['content'];size=item['size'];parts=4
 def fetch(i):
  lo=size*i//parts;hi=size*(i+1)//parts;path=folder/f'{i}.bin';start=path.stat().st_size if path.exists() else 0;assert start<=hi-lo
  if start<hi-lo:
   with requests.get(url,headers={'Range':f'bytes={lo+start}-{hi-1}'},stream=True,timeout=(30,60)) as response:
    response.raise_for_status();assert response.status_code==206 and response.headers['Content-Range']==f'bytes {lo+start}-{hi-1}/{size}'
    with path.open('ab' if start else 'wb') as f:
     done=start;last=start
     for block in response.iter_content(4*1024**2):
      if not block:continue
      f.write(block);done+=len(block)
      if done-last>=128*1024**2:print('CD4 D2 part',i,done,'/',hi-lo,flush=True);last=done
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
 (A/'cd4_D2_download.json').write_text(json.dumps({'file':item,'bytes':size,'md5_verified':True,'metadata_only_selection':'D2 Rest has 15 competition targets>=400 cells, versus D1 Rest2, D1 Stim48hr2, D4 Rest1; no treatment expression score used.','download_script_sha256':b.sha256(Path(__file__))},indent=2),encoding='utf-8');print('CD4 D2 raw source MD5 verified',flush=True)
if __name__=='__main__':main()
