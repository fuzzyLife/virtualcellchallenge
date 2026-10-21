"""Verify staged diagnostic-release provenance and preservation of all earlier protocols."""
from pathlib import Path
import ast,json,hashlib,subprocess
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v90_interval_diagnostic/analysis'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def git(*args):return subprocess.run(['git','-c','safe.directory='+ROOT.as_posix(),*args],cwd=ROOT,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE).stdout
checked=[]
for folder in ['v85_guide_contrast','v86_joint_posterior','v87_conditional_adapter','v88_mixture_quantiles','v89_v86_projection','v90_interval_diagnostic']:
 plan=read(ROOT/'experiments'/folder/'analysis/plan.json')
 for name,h in plan['code_sha256'].items():
  assert sha(ROOT/name)==h,name;assert hashlib.sha256(git('show',':'+name)).hexdigest()==h,name;checked.append(name)
for name in ['.gitattributes','.gitignore']:assert (ROOT/name).read_bytes().startswith(git('show','HEAD:'+name)),name
assert (ROOT/'README.md').read_bytes().endswith(git('show','HEAD:README.md'))
files=[n for n in git('diff','--cached','--name-only','-z').decode().split('\0') if n]
for name in files:
 assert Path(name).suffix.lower() not in {'.npz','.npy','.pt','.vcc','.h5ad','.zip','.bin','.parquet','.pdf','.log'},name
 assert not name.startswith(('bolero/','pie/','tahoe-x1-upstream/')),name
 if '/' not in name and name.endswith('.py'):ast.parse((ROOT/name).read_text(encoding='utf-8-sig'),filename=name)
manifest=read(ROOT/'experiments/v69_count_bridge/analysis/release_manifest.json');text=0
for r in manifest['files']:
 p=ROOT/r['path'];assert p.stat().st_size==r['bytes'],r['path']
 if p.suffix.lower() in ['.py','.json','.md','.csv','.txt']:assert sha(p)==r['sha256'],r['path'];text+=1
for folder in ['v88_mixture_quantiles','v89_v86_projection']:
 assert not read(ROOT/'experiments'/folder/'analysis/release_decision.json')['experimental_release_allowed'];assert not list((ROOT/'experiments'/folder).glob('prediction.vcc'))
plan=read(A/'plan.json');assert not plan['strong_support'];assert plan['status']=='POST_HOC_DIAGNOSTIC_NOT_VALIDATED_IMPROVEMENT';assert read(A/'release_decision.json')['prior_V88_V89_gates_failed']
git('-c','core.whitespace=cr-at-eol','diff','--cached','--check')
ready=(A/'independent_validation.json').exists() and read(A/'independent_validation.json')['passed']
result=dict(passed=True,frozen_source_paths=sorted(set(checked)),README_history_preserved=True,ignore_attribute_history_preserved=True,large_arrays_and_clones_not_staged=True,root_source_syntax_passed=True,V69_inventory_sizes_checked=len(manifest['files']),V69_text_code_hashes_checked=text,V69_large_binary_hashes_rechecked=False,failed_original_V88_V89_decisions_preserved=True,V90_post_hoc_diagnostic_status_explicit=True,upload_integrity_validation_completed=bool(ready),staged_paths_before_receipt=len(files),auditor_sha256=sha(__file__))
(A/'repository_audit.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf-8'));print('V90 repository/preservation audit passed',len(set(checked)),'frozen source paths; upload validated',bool(ready),flush=True)
