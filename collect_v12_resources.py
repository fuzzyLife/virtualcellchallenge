"""Collect pinned public Flex source bins and an isolated official scorer checkout."""
from pathlib import Path
import argparse, hashlib, json, time, zipfile, io
from concurrent.futures import ThreadPoolExecutor
import requests
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'reference/v12_resources'
PIN='0ce25cf4c495e76e60dfca4ddc8f822507aa0e18'

def digest(path,algorithm='sha256'):
 h=hashlib.new(algorithm)
 with path.open('rb') as f:
  for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
 return h.hexdigest()

def download(record):
 path=OUT/'viperturb'/record['key'];path.parent.mkdir(parents=True,exist_ok=True)
 expected=record['checksum'].split(':')[1]
 if path.exists():
  assert path.stat().st_size==record['size'] and digest(path,'md5')==expected
  print('Verified',path.name,flush=True);return {'file':str(path.relative_to(ROOT)),'md5':expected,'sha256':digest(path),'bytes':path.stat().st_size,'url':record['links']['self']}
 partial=path.with_suffix(path.suffix+'.partial')
 for attempt in range(5):
  offset=partial.stat().st_size if partial.exists() else 0
  if offset==record['size']:break
  response=requests.get(record['links']['self'],headers={'Range':f'bytes={offset}-'} if offset else {},stream=True,timeout=(60,120))
  response.raise_for_status()
  if offset and response.status_code!=206:offset=0
  mode='ab' if offset else 'wb';start=time.time();last=start;written=offset
  try:
   with partial.open(mode) as f:
    for block in response.iter_content(4*1024**2):
     if not block:continue
     f.write(block);written+=len(block)
     if time.time()-last>30:
      print(path.name,f'{written/2**30:.2f}/{record["size"]/2**30:.2f} GiB',f'{(written-offset)/2**20/max(time.time()-start,1):.1f} MiB/s',flush=True);last=time.time()
   break
  except (requests.RequestException,OSError) as error:
   print('Retry',path.name,attempt,type(error).__name__,flush=True)
 if partial.stat().st_size!=record['size']:raise ValueError('Incomplete '+path.name)
 assert digest(partial,'md5')==expected,path.name
 partial.replace(path)
 result={'file':str(path.relative_to(ROOT)),'md5':expected,'sha256':digest(path),'bytes':path.stat().st_size,'url':record['links']['self']}
 side=path.with_suffix(path.suffix+'.json');side.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 print('Downloaded and checksum-verified',path.name,flush=True);return result

def scorer():
 path=OUT/'scorer';url=f'https://codeload.github.com/ArcInstitute/cell-eval2/zip/{PIN}'
 if not (path/'pyproject.toml').exists():
  r=requests.get(url,timeout=120);r.raise_for_status()
  with zipfile.ZipFile(io.BytesIO(r.content)) as archive:
   prefix=archive.namelist()[0].split('/')[0]+'/'
   for item in archive.infolist():
    if item.is_dir():continue
    name=item.filename.removeprefix(prefix)
    if name.startswith(('src/','docs/','tests/')) or name in ['LICENSE','README.md','pyproject.toml']:
     target=path/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(archive.read(item))
 old=json.loads((ROOT/'analysis/v8_official_metric_review.json').read_text())
 for name,record in old['files'].items():assert digest(path/name)==record['sha256'],name
 print('Official scorer ready',PIN,flush=True)

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--sources',action='store_true');args=parser.parse_args()
 OUT.mkdir(parents=True,exist_ok=True);scorer()
 if args.sources:
  record=json.loads((ROOT/'analysis/post_v11_research/data/viperturb_record.json').read_text())
  files=[f for f in record['files'] if f['key'] in ['genome_wide_binA.RDS','genome_wide_binB.RDS','genome_wide_binC.RDS']]
  with ThreadPoolExecutor(max_workers=3) as pool:results=list(pool.map(download,files))
  (OUT/'viperturb_manifest.json').write_text(json.dumps({'doi':record['doi'],'license':record['metadata']['license'],'files':results},indent=2)+'\n',encoding='utf-8')
