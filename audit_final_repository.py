"""Verify final release provenance and preserve historical scientific protocols/readiness records."""
from pathlib import Path
import hashlib,json,subprocess,ast
R=Path(__file__).resolve().parent;A=R/'experiments/final_v91/analysis'
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def git(*args):return subprocess.run(['git','-c','safe.directory='+R.as_posix(),*args],cwd=R,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE).stdout
def main():
 paths=set()
 for folder in ['v85_guide_contrast','v86_joint_posterior','v87_conditional_adapter','v88_mixture_quantiles','v89_v86_projection','v90_interval_diagnostic','final_v91']:
  plan=read(R/'experiments'/folder/'analysis/plan.json')
  for name,h in plan['code_sha256'].items():assert sha(R/name)==h,name;assert hashlib.sha256(git('show',':'+name)).hexdigest()==h,name;paths.add(name)
 plan=read(A/'plan.json')
 for name,h in plan['source_receipt_sha256'].items():assert sha(R/name)==h,name
 inputs=read(A/'input_audit.json')['inputs']
 for name,row in inputs.items():assert (R/name).stat().st_size==row['bytes'] and sha(R/name)==row['sha256'],name
 baseline='dff10f1'
 assert (R/'README.md').read_bytes().endswith(git('show',baseline+':README.md'))
 for name in ['.gitignore','.gitattributes']:assert (R/name).read_bytes().startswith(git('show',baseline+':'+name))
 manifest=read(R/'experiments/v69_count_bridge/analysis/release_manifest.json');texts=0
 for row in manifest['files']:
  p=R/row['path'];assert p.stat().st_size==row['bytes'],row['path']
  if p.suffix.lower() in ['.py','.json','.md','.csv','.txt']:assert sha(p)==row['sha256'],row['path'];texts+=1
 for folder in ['v88_mixture_quantiles','v89_v86_projection']:assert not read(R/'experiments'/folder/'analysis/release_decision.json')['experimental_release_allowed']
 ready=read(A/'submission_ready.json');ind=read(A/'independent_validation.json');assert ind['passed'];assert sha(R/'experiments/final_v91/prediction.vcc')==ready['sha256']==ind['submission_sha256']
 staged=[s for s in git('diff','--cached','--name-only','-z').decode().split('\0') if s]
 for name in staged:
  assert Path(name).suffix.lower() not in ['.npz','.npy','.pt','.vcc','.h5ad','.zip','.bin','.parquet','.pdf','.log'],name
  assert not name.startswith(('bolero/','pie/','tahoe-x1-upstream/')),name
  if '/' not in name and name.endswith('.py'):ast.parse((R/name).read_text(encoding='utf-8-sig'))
 git('-c','core.whitespace=cr-at-eol','diff','--cached','--check')
 result=dict(passed=True,frozen_code_files=len(paths),code_file_and_git_index_hashes_match=True,all_final_input_files_unchanged=len(inputs),all_source_receipt_hashes_match=True,README_history_preserved=True,old_ignore_attributes_preserved=True,V69_inventory_sizes_checked=len(manifest['files']),V69_text_code_hashes_checked=texts,V69_large_binary_hashes_rechecked=False,V88_V89_failures_preserved=True,final_upload_sha256=ready['sha256'],official_prep_and_independent_validation_passed=True,large_arrays_and_unrelated_clones_not_staged=True,auditor_sha256=sha(Path(__file__)))
 (A/'repository_audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print('Final repository and preservation audit passed',len(paths),'frozen code files',flush=True)
if __name__=='__main__':main()
