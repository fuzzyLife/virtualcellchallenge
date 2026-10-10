"""Verified, bounded HTTP-range collection of the frozen VIPerturb resource."""
import hashlib,json,time,threading,os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import requests
ROOT=Path(__file__).resolve().parent
D=ROOT/'reference/v12_resources/viperturb'
BLOCK=8*1024**2
lock=threading.Lock();finished=0

def main(workers=48):
 global finished
 record=json.loads((ROOT/'analysis/post_v11_research/data/viperturb_record.json').read_text())
 files=[f for f in record['files'] if f['key'] in ['genome_wide_binA.RDS','genome_wide_binB.RDS','genome_wide_binC.RDS']]
 jobs=[];D.mkdir(parents=True,exist_ok=True)
 for f in files:
  for start in range(0,f['size'],BLOCK):jobs.append((f,start,min(start+BLOCK,f['size'])-1))
 def fetch(job):
  global finished
  f,start,end=job;folder=D/(f['key']+'.ranges');folder.mkdir(exist_ok=True)
  out=folder/f'{start:012d}.bin'
  if not (out.exists() and out.stat().st_size==end-start+1):
   for attempt in range(8):
    try:
     r=requests.get(f['links']['self'],headers={'Range':f'bytes={start}-{end}'},stream=True,timeout=(30,180))
     if r.status_code in [429,503]:
      time.sleep(min(120,float(r.headers.get('Retry-After',60))));continue
     r.raise_for_status()
     assert r.status_code==206 and r.headers['Content-Range']==f"bytes {start}-{end}/{f['size']}"
     with out.with_suffix('.partial').open('wb') as h:
      for chunk in r.iter_content(1024**2):h.write(chunk)
     assert out.with_suffix('.partial').stat().st_size==end-start+1
     os.replace(out.with_suffix('.partial'),out);break
    except Exception:
     if attempt==7:raise
     time.sleep(2+attempt*3)
   else:raise RuntimeError('HTTP range retries exhausted')
  with lock:
   finished+=1
   if finished%12==0:print(f'Verified ranges {finished}/{len(jobs)}',flush=True)
 with ThreadPoolExecutor(max_workers=workers) as ex:
  for task in as_completed([ex.submit(fetch,j) for j in jobs]):task.result()
 manifest=[]
 for f in files:
  dest=D/f['key'];md5=hashlib.md5();sha=hashlib.sha256()
  with dest.open('wb') as h:
   for start in range(0,f['size'],BLOCK):
    chunk=(D/(f['key']+'.ranges')/f'{start:012d}.bin').read_bytes();h.write(chunk);md5.update(chunk);sha.update(chunk)
  assert dest.stat().st_size==f['size'] and 'md5:'+md5.hexdigest()==f['checksum']
  manifest.append({'file':f['key'],'bytes':f['size'],'md5':md5.hexdigest(),'sha256':sha.hexdigest()})
 (D.parent/'viperturb_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print('All bins verified',flush=True)
if __name__=='__main__':main()
