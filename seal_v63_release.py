"""Seal V62/V63 implementation, results and artifact receipts without uploading."""
from pathlib import Path
import argparse
from v57_pie import read,write,sha
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v63_identity_guard/analysis'

def main(check=False):
 path=A/'release_manifest.json'
 if check:
  r=read(path)
  for n,h in r['files'].items():assert sha(ROOT/n)==h,n
  print('V62/V63 release seal verified',len(r['files']),'files',flush=True);return
 assert read(A/'independent_validation.json')['passed'];selected=set()
 for pattern in ('*v62*.py','*v63*.py','*v64*.py','V62*.md','V63*.md','V64*.md'):selected.update(ROOT.glob(pattern))
 for folder in ('v62_conserved_response','v63_identity_guard','v64_response_magnitude'):
  for p in (ROOT/'experiments'/folder/'analysis').rglob('*'):
   if p.is_file() and p.suffix in ('.json','.csv','.md','.txt','.png') and p.name not in ('release_manifest.json','prediction_control_rows.csv'):selected.add(p)
 write(path,dict(files={str(p.relative_to(ROOT)):sha(p) for p in sorted(selected)},uploaded=False,scope='V62/V63/V64 root methods/protocol/reports and compact scientific receipts. Numerical arrays and large package/count files are bound by their explicit receipt hashes. README is a mutable status overview and excluded.'))
 print('Sealed V62/V63',len(selected),'files',flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');main(p.parse_args().check)
