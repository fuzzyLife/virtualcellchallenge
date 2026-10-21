"""Pinned experimental-response source acquisition for V83; no model inference."""
from pathlib import Path
import argparse,hashlib,json,requests
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v83_context_expansion/analysis';R=ROOT/'reference/v83_sources'
def main():
 plan=json.loads((A/'PIE_replogle_nadig_essential_inventory.json').read_text());rev=plan['revision'];R.mkdir(parents=True,exist_ok=True);rows=[f for f in plan['files'] if f['path'] in ['de/hepg2.parquet','de/jurkat.parquet','de/RPE1.parquet','de/K562.parquet']]
 def fetch(f):
  p=R/Path(f['path']).name;expected=f.get('lfs',{}).get('oid');url=f'https://huggingface.co/datasets/arcinstitute/PIE_replogle_nadig_essential/resolve/{rev}/{f["path"]}'
  if not p.exists():
   pending=p.with_suffix(p.suffix+'.partial');assert not pending.exists();h=hashlib.sha256()
   with requests.get(url,stream=True,timeout=(20,90)) as response:
    response.raise_for_status()
    with pending.open('xb') as stream:
     for block in response.iter_content(4*1024**2):
      if block:stream.write(block);h.update(block)
   assert pending.stat().st_size==f['size'];actual=h.hexdigest();assert not expected or expected==actual,(expected,actual);pending.replace(p)
  else:
   with p.open('rb') as stream:actual=hashlib.file_digest(stream,'sha256').hexdigest()
   assert p.stat().st_size==f['size'];assert not expected or actual==expected
  receipt=dict(repository='arcinstitute/PIE_replogle_nadig_essential',revision=rev,path=f['path'],bytes=p.stat().st_size,sha256=actual,expected_lfs_sha256=expected,publisher_digest_verified=bool(expected),url=url,collector_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());p.with_suffix('.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8');print('Downloaded measured DE',f['path'],p.stat().st_size,flush=True);return receipt
 with ThreadPoolExecutor(max_workers=3) as pool:records=list(pool.map(fetch,rows))
 (A/'source_downloads.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
if __name__=='__main__':main()
