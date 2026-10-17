"""Collect official LINCS genetic-consensus metadata before choosing response features."""
import gzip,hashlib,json
from pathlib import Path
import requests
import build_submission as b
ROOT=b.ROOT;OUT=ROOT/'experiments/v43_lincs_response';A=OUT/'analysis';R=OUT/'reference'
BASE='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106127/suppl/'
FILES=['GSE106127_SHA512SUMS.txt.gz','GSE106127_CGS_PRIME_meta.txt.gz','GSE106127_gene_info.txt.gz','GSE106127_README.pdf.gz','GSE106127_CGS_PRIME_n33839x978.gctx.gz']
def fetch(name):
 R.mkdir(parents=True,exist_ok=True);p=R/name
 if not p.exists():
  with requests.get(BASE+name,stream=True,timeout=(30,90)) as r:
   r.raise_for_status();total=int(r.headers.get('Content-Length','0'));assert total<500000000 if total else True
   pending=p.with_suffix(p.suffix+'.partial');n=0
   with pending.open('wb') as f:
    for block in r.iter_content(1048576):n+=len(block);assert n<500000000;f.write(block)
  if total:assert total==n
  pending.replace(p)
 return p
def main():
 A.mkdir(parents=True,exist_ok=True);fetch(FILES[0]);checks=gzip.decompress((R/FILES[0]).read_bytes()).decode();(A/'publisher_SHA512SUMS.txt').write_text(checks,encoding='utf-8');records=[]
 for name in FILES[1:]:
  p=fetch(name);payload=p.read_bytes();h=hashlib.sha512(payload).hexdigest()
  matches=[line for line in checks.splitlines() if name in line]
  assert len(matches)==1,(name,matches);assert h==matches[0].split()[0],name
  records.append({'file':name,'url':BASE+name,'bytes':len(payload),'sha512':h,'publisher_sha512_verified':True,'sha256':hashlib.sha256(payload).hexdigest()});print('LINCS verified',name,len(payload),flush=True)
 b.write_json(A/'collection.json',{'files':records,'script_sha256':b.sha256(Path(__file__)),'accession':'GSE106127','scope':'CGS PRIME genetic consensus landmark profiles, source units to be inspected before any fit. No count or fold-change conversion implied.','publisher_checksum_manifest_sha256':b.sha256(R/FILES[0])})
 import pandas as pd
 meta=pd.read_csv(R/FILES[1],sep='\t');genes=pd.read_csv(R/FILES[2],sep='\t');_,axis,targets=b.read_spec(ROOT)
 print('metadata shape',meta.shape,'columns',meta.columns.tolist(),flush=True);print(meta.head(3).to_string(index=False),flush=True);print('genes',genes.head(3).to_string(index=False),flush=True)
 for col in meta.columns:
  if col.lower() in ['pert_iname','pert_id','gene','gene_name','pert_type','cell_id']:print(col,'target overlap',len(set(meta[col].astype(str))&set(targets)),'unique',meta[col].nunique(),flush=True)
if __name__=='__main__':main()

