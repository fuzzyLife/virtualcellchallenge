"""Idempotently unpack checksum-verified public LINCS payloads."""
import gzip,json,shutil
from pathlib import Path
import build_submission as b
R=b.ROOT/'experiments/v43_lincs_response/reference';A=R.parent/'analysis'
def main():
 rec=json.loads((A/'collection.json').read_text())
 output=[]
 for item in rec['files']:
  if not item['file'].endswith(('.gctx.gz','.pdf.gz')):continue
  src=R/item['file'];assert b.sha256(src)==item['sha256'];dst=src.with_suffix('')
  # Independently verify existing decompression by streamed uncompressed SHA.
  import hashlib
  h=hashlib.sha256()
  with gzip.open(src,'rb') as f:
   while block:=f.read(8*1024*1024):h.update(block)
  if not dst.exists():
   with gzip.open(src,'rb') as f,open(dst,'wb') as o:shutil.copyfileobj(f,o,8*1024*1024)
  assert b.sha256(dst)==h.hexdigest()
  output.append({'file':dst.name,'sha256':h.hexdigest(),'bytes':dst.stat().st_size})
 b.write_json(A/'unpacked.json',{'files':output,'script_sha256':b.sha256(Path(__file__))})
if __name__=='__main__':main()

