"""Verify or fetch the small pinned research inputs; never execute downloaded code."""
from pathlib import Path
import argparse,json,hashlib
import requests
import build_submission as b
A=b.ROOT/'analysis/post_v17_research'
def main(fetch=False):
 plan=json.loads((A/'source_manifest.json').read_text());verified=[]
 for item in plan['sources']:
  path=A/item['file'];assert path.resolve().parent==A.resolve()
  if not path.exists():
   if not fetch:raise FileNotFoundError(path)
   if item['kind']=='api_snapshot':raise ValueError('Historical API snapshot must be restored from Git: '+item['file'])
   r=requests.get(item['url'],timeout=120);r.raise_for_status();assert len(r.content)==item['bytes'] and hashlib.sha256(r.content).hexdigest()==item['sha256'];path.write_bytes(r.content)
  assert path.stat().st_size==item['bytes'] and b.sha256(path)==item['sha256'],str(path)
  verified.append(item['file'])
 b.write_json(A/'source_verification.json',{'passed':True,'verified_files':len(verified),'manifest_sha256':b.sha256(A/'source_manifest.json'),'collector_sha256':b.sha256(Path(__file__)),'downloaded_author_code_executed':False,'model_weights_downloaded':False});print('Pinned research source checks passed:',len(verified),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--fetch-missing',action='store_true');main(p.parse_args().fetch_missing)
