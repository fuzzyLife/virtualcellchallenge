"""Collect pinned ungated TabICLv2 regression weights, with verified byte ranges."""
import requests,hashlib,json
from pathlib import Path
from urllib.parse import urljoin
ROOT=Path(__file__).resolve().parent;D=ROOT/'reference/v33_tabular';D.mkdir(parents=True,exist_ok=True)
REV='4dcd344ece2c00be9e831fdd35bed57b5ad83e19';NAME='tabicl-regressor-v2-20260212.ckpt';SHA='0db9cb538f114e79026bf08f45f41ad8dd7ad2de2aaca9a5ca8cd3bd9748ae7a';SIZE=114324594
p=D/NAME;url=f'https://huggingface.co/jingang/TabICL/resolve/{REV}/{NAME}';tmp=p.with_suffix('.partial')
if not p.exists():
 start=tmp.stat().st_size if tmp.exists() else 0
 with tmp.open('ab') as f:
  for lo in range(start,SIZE,4*1024*1024):
   hi=min(SIZE-1,lo+4*1024*1024-1);u=url+f'?download=true&v33range={lo}'
   for attempt in range(3):
    try:
     for _ in range(5):
      r=requests.get(u,headers={'Range':f'bytes={lo}-{hi}'},stream=True,allow_redirects=False,timeout=(20,30))
      if r.status_code in [301,302,303,307,308]:u=urljoin(u,r.headers['Location']);r.close();continue
      break
     r.raise_for_status();assert r.status_code==206 and r.headers['Content-Range']==f'bytes {lo}-{hi}/{SIZE}';block=b''.join(r.iter_content(1024*1024));r.close();assert len(block)==hi-lo+1;f.write(block);f.flush();print('TabICL verified range',hi+1,'/',SIZE,flush=True);break
    except Exception:
     if attempt==2:raise
     u=url+f'?download=true&v33range={lo}&retry={attempt+1}'
 tmp.replace(p)
assert p.stat().st_size==SIZE;assert hashlib.file_digest(p.open('rb'),'sha256').hexdigest()==SHA
receipt={'url':url,'revision':REV,'bytes':SIZE,'sha256':SHA,'license':'BSD-3-Clause','gated':False,'packages':{'tabicl':'2.2.0','einops':'0.8.2','catboost':'1.2.10'},'weights_loader':'Official TabICL uses torch.load(weights_only=True); no pickled foreign model class.'};(D/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print('Checkpoint SHA256 verified',SHA,flush=True)