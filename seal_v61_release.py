"""Seal scientific code and compact state-release receipts; verify without mutation."""
import argparse
from v57_pie import ROOT,sha,read,write
A=ROOT/'experiments/v61_pooled_state_response/analysis'

def files():
 selected=set()
 for pattern in ['*v59*.py','*v60*.py','*v61*.py','V59*.md','V60*.md','V61*.md']:
  selected.update(p for p in ROOT.glob(pattern) if p.is_file())
 selected.update(ROOT/name for name in ['document_state_release_review.py','evaluate_state_population_diagnostic.py'])
 for folder in ['v59_state_composition','v60_depth_adjusted_states','v61_pooled_state_response']:
  base=ROOT/'experiments'/folder/'analysis'
  for path in base.rglob('*'):
   if path.is_file() and path.suffix in ['.json','.csv','.png'] and path.name not in ['release_manifest.json','prediction_control_rows.csv'] and 'arrays' not in path.relative_to(base).parts:selected.add(path)
 return sorted(selected)

def main(check=False):
 target=A/'release_manifest.json'
 if check:
  p=read(target)
  for name,h in p['files'].items():assert sha(ROOT/name)==h,name
  print('State release seal verified',len(p['files']),'files',flush=True)
 else:
  assert read(A/'independent_validation.json')['passed'];write(target,dict(files={str(p.relative_to(ROOT)):sha(p) for p in files()},seal_sha256=sha(__file__),scope='Root V59-V61 implementation/protocol/reports and compact analysis files. Large numerical arrays are covered by model/export manifests; the final upload package is covered by independent_validation.json.',uploaded=False));print('State release sealed',len(read(target)['files']),'files',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');main(p.parse_args().check)
