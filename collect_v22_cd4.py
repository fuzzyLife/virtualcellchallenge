"""Fetch one author-published donor-resolved CD4 count file with MD5 checks."""
from pathlib import Path
import argparse,hashlib,json,time,requests
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v22_state_response/analysis';R=ROOT/'reference/v22_cd4'
def main():
 p=argparse.ArgumentParser();p.add_argument('--download',action='store_true');args=p.parse_args();A.mkdir(parents=True,exist_ok=True);R.mkdir(parents=True,exist_ok=True)
 u='https://zenodo.org/api/records/21729034';saved=A/'cipher_zenodo_record.json'
 if not saved.exists():
  r=requests.get(u,timeout=40);r.raise_for_status();saved.write_text(json.dumps(r.json(),indent=2),encoding='utf-8')
 record=json.loads(saved.read_text(encoding='utf-8'));item=next(v for v in record['files'] if v['key']=='Marson2025_D1_Rest_filtered.h5ad');dest=R/item['key'];expected=item['checksum'].split(':')[1]
 provenance={'api_url':u,'record_id':record['id'],'file':item,'metadata_sha256':hashlib.sha256(saved.read_bytes()).hexdigest(),'role':'Candidate donor-resolved CD4 raw-cell benchmark; same Marson study may already supply source summaries. No new independent study claim.'}
 if not args.download:print(json.dumps(provenance,indent=2));return
 if not dest.exists():
  part=dest.with_suffix('.h5ad.partial');start=part.stat().st_size if part.exists() else 0;url=item['links'].get('self') or item['links']['content'];headers={'Range':f'bytes={start}-'} if start else {}
  with requests.get(url,headers=headers,stream=True,timeout=(30,60)) as response:
   response.raise_for_status()
   if start and response.status_code!=206:raise RuntimeError('Server refused safe resume')
   with part.open('ab' if start else 'wb') as f:
    done=start;last=done;t=time.monotonic()
    for chunk in response.iter_content(4*1024*1024):
     if not chunk:continue
     f.write(chunk);done+=len(chunk)
     if done-last>=64*1024*1024:print('CD4 download',done,'/',item['size'],'bytes',round((done-start)/max(time.monotonic()-t,1)/1e6,2),'MB/s',flush=True);last=done
  assert part.stat().st_size==item['size']
  with part.open('rb') as f:assert hashlib.file_digest(f,'md5').hexdigest()==expected
  part.replace(dest)
 with dest.open('rb') as f:assert hashlib.file_digest(f,'md5').hexdigest()==expected
 provenance.update(downloaded=True,bytes=dest.stat().st_size,md5_verified=True)
 (A/'cd4_download.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8');print('CD4 original count file verified',flush=True)
if __name__=='__main__':main()
