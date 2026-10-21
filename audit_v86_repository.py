"""Check staged scientific bytes, retained README history and excluded large outputs."""
from pathlib import Path
import ast,hashlib,json,subprocess
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v86_joint_posterior/analysis'
def git(*args):
 return subprocess.run(['git','-c','safe.directory='+ROOT.as_posix(),*args],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=ROOT).stdout
checked=[]
for version,folder in [(85,'v85_guide_contrast'),(86,'v86_joint_posterior')]:
 plan=json.loads((ROOT/f'experiments/{folder}/analysis/plan.json').read_text())
 for name,expected in plan['code_sha256'].items():
  actual=hashlib.sha256(git('show',':'+name)).hexdigest();assert actual==expected,name;checked.append(name)
old=git('show','HEAD:README.md');assert (ROOT/'README.md').read_bytes().endswith(old),'Prior README history changed'
files=[f for f in git('diff','--cached','--name-only','-z').decode().split('\0') if f]
for name in files:
 assert Path(name).suffix.lower() not in {'.npz','.npy','.h5ad','.vcc','.log','.parquet','.bin','.zip'},name
 assert not ('v85_guide_contrast/analysis/teachers/' in name or 'v85_guide_contrast/analysis/linear_teachers/' in name),name
 if '/' not in name and name.endswith('.py'):ast.parse((ROOT/name).read_text(encoding='utf-8-sig'),filename=name)
git('-c','core.whitespace=cr-at-eol','diff','--cached','--check')
result=dict(passed=True,frozen_code_paths=sorted(set(checked)),staged_paths=len(files),large_generated_arrays_and_repeated_full_traces_excluded=True,root_script_syntax_checked=True,README_history_bytes_preserved=True,whitespace_check_passed=True,auditor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(A/'repository_audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print('V85/V86 repository evidence audit passed',len(set(checked)),'frozen code paths',flush=True)
