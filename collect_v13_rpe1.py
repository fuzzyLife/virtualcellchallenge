"""Download the author's raw RPE1 pseudobulk release; verify its published MD5."""
from pathlib import Path
import hashlib,json,requests
import build_submission as build
ROOT=build.ROOT
D=ROOT/'reference/v13_resources';D.mkdir(parents=True,exist_ok=True)
URL='https://ndownloader.figshare.com/files/35775581'
EXPECTED='74765fa87635467a869ea972356ae0e7'
SIZE=95350546

def main():
 path=D/'rpe1_raw_bulk_01.h5ad'
 if not path.exists():
  with requests.get(URL,stream=True,timeout=90) as response:
   response.raise_for_status()
   with path.with_suffix('.partial').open('wb') as f:
    for block in response.iter_content(8*1024*1024):f.write(block)
  temporary=path.with_suffix('.partial');assert temporary.stat().st_size==SIZE
  assert hashlib.md5(temporary.read_bytes()).hexdigest()==EXPECTED
  temporary.replace(path)
 assert path.stat().st_size==SIZE and hashlib.md5(path.read_bytes()).hexdigest()==EXPECTED
 build.write_json(D/'rpe1_download.json',{'url':URL,'bytes':SIZE,'md5':EXPECTED,'sha256':build.sha256(path),
  'dataset':'Replogle et al. 2022, RPE1 essential CRISPRi, raw pseudobulk',
  'record':'https://plus.figshare.com/articles/dataset/20029387','license':'CC BY 4.0','script_sha256':build.sha256(Path(__file__))})
 import anndata as ad
 data=ad.read_h5ad(path)
 print(data,flush=True);print(data.obs.head().to_string(),flush=True);print(data.var.head().to_string(),flush=True);print('layers',list(data.layers),'uns',list(data.uns),flush=True)
if __name__=='__main__':main()
