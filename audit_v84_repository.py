"""Check that Git preserves frozen scientific bytes and excludes large outputs."""
from pathlib import Path
import hashlib,json,subprocess,ast
ROOT=Path(__file__).resolve().parent
A=ROOT/'experiments/v84_intervention_operator/analysis'
SAFE='safe.directory='+ROOT.as_posix()
def git(*args):
 return subprocess.run(['git','-c',SAFE,*args],cwd=ROOT,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE).stdout
protocol=json.loads((A/'execution_protocol.json').read_text())
checked=[]
for name,expected in protocol['code_sha256'].items():
 raw=git('show',':'+name);assert hashlib.sha256(raw).hexdigest()==expected,name
 checked.append(name)
for path in (ROOT/'experiments/v83_context_expansion/analysis/author_code').iterdir():
 if path.is_file():
  name=path.relative_to(ROOT).as_posix();assert git('show',':'+name)==path.read_bytes(),name;checked.append(name)
old=git('show','HEAD:README.md');assert (ROOT/'README.md').read_bytes().endswith(old),'README history was rewritten'
files=git('diff','--cached','--name-only','-z').decode().split('\0');files=[f for f in files if f]
for name in files:
 assert Path(name).suffix.lower() not in {'.npz','.npy','.h5ad','.vcc','.log','.bin','.parquet','.zip'},name
 if '/' not in name and name.endswith('.py'):ast.parse((ROOT/name).read_text(encoding='utf-8-sig'),filename=name)
git('-c','core.whitespace=cr-at-eol','diff','--cached','--check')
result=dict(passed=True,frozen_source_index_bytes_checked=checked,staged_paths=len(files),large_generated_outputs_not_staged=True,README_history_bytes_preserved=True,root_python_syntax_checked=True,whitespace_check_passed=True,frozen_source_line_endings_preserved_with_explicit_attributes=True,auditor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(A/'repository_audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print('V84 repository evidence audit passed',len(checked),'frozen source files;',len(files),'staged paths',flush=True)
