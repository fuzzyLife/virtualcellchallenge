"""Seal compact V65/V66 scientific outputs while preserving all V63 files."""
from pathlib import Path
import argparse
from v57_pie import read,write,sha
ROOT=Path(__file__).resolve().parent;BASE=ROOT/'experiments/v65_closed_response/analysis'

def main(check=False):
 decision=read(BASE/'release_selection.json');version=decision['selected_version'];assert version in (65,66);folder='v65_closed_response' if version==65 else 'v66_assay_strength';a=ROOT/'experiments'/folder/'analysis';path=a/'release_manifest.json';prior=ROOT/'experiments/v63_identity_guard/analysis/release_manifest.json'
 for n,h in read(prior)['files'].items():assert sha(ROOT/n)==h,n
 if check:
  r=read(path);assert r['prior_V63_manifest_sha256']==sha(prior)
  for n,h in r['files'].items():assert sha(ROOT/n)==h,n
  assert (a.parent/'prediction.vcc').stat().st_size==read(a/'submission_ready.json')['bytes']
  print('V65/V66 compact seal verified',len(r['files']),'files; prior V63 scientific files unchanged.',flush=True);return
 assert not path.exists();assert read(a/'independent_validation.json')['passed'] and read(a/'count_replay_audit.json')['passed'];chosen=set()
 for pattern in ('*v65*.py','*v66*.py','V65*.md','V66*.md'):chosen.update(ROOT.glob(pattern))
 for dirname in ('v65_closed_response','v66_assay_strength'):
  for p in (ROOT/'experiments'/dirname/'analysis').rglob('*'):
   if p.is_file() and p.suffix in ('.json','.csv','.md','.txt','.png') and p.name not in ('release_manifest.json','prediction_control_rows.csv'):chosen.add(p)
 write(path,dict(files={str(p.relative_to(ROOT)):sha(p) for p in sorted(chosen)},selected_version=version,prior_V63_manifest_sha256=sha(prior),package_sha256=read(a/'submission_ready.json')['sha256'],uploaded=False,private_score=None,scope='New scientific implementation/protocol/reports and compact receipts, including negative alternatives. Large model/data/count/package arrays are bound through their receipt hashes. README is mutable and excluded. Future publication receipts must be added separately without modifying this pre-upload seal.'))
 print('Sealed V65/V66',len(chosen),'compact scientific files; selected V'+str(version),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');main(p.parse_args().check)
