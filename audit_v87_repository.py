"""Audit staged V87 source bytes and earlier scientific records without rewriting old receipts."""
from pathlib import Path
import ast,hashlib,json,subprocess
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v87_conditional_adapter/analysis'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def git(*args):return subprocess.run(['git','-c','safe.directory='+ROOT.as_posix(),*args],cwd=ROOT,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE).stdout
checked=[]
for folder in ['v85_guide_contrast','v86_joint_posterior','v87_conditional_adapter']:
 for name,h in read(ROOT/'experiments'/folder/'analysis/plan.json')['code_sha256'].items():
  assert sha(ROOT/name)==h,name;assert hashlib.sha256(git('show',':'+name)).hexdigest()==h,name;checked.append(name)
assert (ROOT/'README.md').read_bytes().endswith(git('show','HEAD:README.md'))
for name in ['.gitattributes','.gitignore']:assert (ROOT/name).read_bytes().startswith(git('show','HEAD:'+name)),name
files=[f for f in git('diff','--cached','--name-only','-z').decode().split('\0') if f]
explicit={'.gitattributes','.gitignore','README.md','BOLERO_REVIEW_V87.md','V87_RESULTS.md'}
for name in files:
 assert name in explicit or (('/' not in name) and ('v87' in name) and name.endswith('.py')) or name.startswith('experiments/v87_conditional_adapter/analysis/'),name
 assert Path(name).suffix.lower() not in {'.npz','.npy','.pt','.h5ad','.vcc','.log','.parquet','.bin','.zip','.pdf'},name
 if '/' not in name and name.endswith('.py'):ast.parse((ROOT/name).read_text(encoding='utf-8-sig'),filename=name)
# V69's old seal inventories all text/code and binary sizes; avoid rehashing multi-GB predictions.
manifest=read(ROOT/'experiments/v69_count_bridge/analysis/release_manifest.json');text_count=0
for r in manifest['files']:
 p=ROOT/r['path'];assert p.stat().st_size==r['bytes'],r['path']
 if p.suffix.lower() in ['.py','.json','.md','.csv','.txt']:assert sha(p)==r['sha256'],r['path'];text_count+=1
for name in ['synthetic_tests.json','model_audit.json','prediction_audit.json','metric_audit.json']:assert read(A/name)['passed'],name
assert not read(A/'mean_results.json')['passes'];decision=read(A/'release_decision.json');assert not decision['new_upload_generated'] and not decision['count_test_run'];assert not list(A.parent.glob('*.vcc'))
git('-c','core.whitespace=cr-at-eol','diff','--cached','--check')
result=dict(passed=True,frozen_code_paths=sorted(set(checked)),README_history_preserved=True,ignore_attribute_history_preserved=True,staged_paths_before_this_receipt=len(files),root_source_syntax_checked=True,V69_old_inventory_sizes_checked=len(manifest['files']),V69_code_text_hashes_checked=text_count,V69_large_binary_hashes_rechecked=False,all_numerical_and_metric_audits_passed=True,failed_candidate_correctly_has_no_count_or_upload=True,only_authorized_V87_paths_staged=True,large_arrays_models_cloned_repo_and_PDF_not_staged=True,auditor_sha256=sha(__file__))
(A/'repository_audit.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf-8'));print('V87 repository/preservation audit passed',len(set(checked)),'frozen source paths',flush=True)
