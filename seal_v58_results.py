"""Seal/check the completed V58 scientific reports and result receipts."""
import argparse
from v58_pie import *

def main(check=False):
 path=A/'result_manifest.json'
 if check:
  manifest=read(path)
  for name,digest in manifest['files'].items():assert sha(ROOT/name)==digest,name
  print('V58 result seal verified',len(manifest['files']),'files',flush=True);return
 assert read(A/'final_audit.json')['passed'] and read(A/'truth_population_replay_audit.json')['passed'] and read(A/'array_exports.json')['passed']
 files=[p for p in A.rglob('*') if p.is_file() and p.suffix in ['.json','.csv','.png'] and p!=path and 'arrays' not in p.relative_to(A).parts]
 files+=list(ROOT.glob('*v58*.py'))+list(ROOT.glob('V58*.md'))
 write(path,{'passed':True,'plan_sha256':sha(A/'plan.json'),'files':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in sorted(set(files))},'large_array_exports_manifest_sha256':sha(A/'array_exports.json'),'new_submission_created':False,'reason':'Candidate failed held-CD4 direction/significance checks; keep V24 champion.','private_competition_score':None})
 print('Sealed V58 results',len(files),'files',flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');main(p.parse_args().check)
