"""Finalize documentation and exact-byte seals after all package checks succeed."""
from pathlib import Path
import subprocess,sys,time,json
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v63_identity_guard/analysis'
while not (A/'independent_validation.json').exists():time.sleep(5)
r=json.loads((A/'independent_validation.json').read_text());assert r['passed']
for script,args in [('document_v63_release.py',[]),('seal_v61_release.py',['--check']),('seal_v63_release.py',[]),('seal_v63_release.py',['--check'])]:
 print('V63 final review',script,*args,flush=True);subprocess.run([sys.executable,'-X','utf8',str(ROOT/script),*args],check=True)
text=(ROOT/'V63_SUBMISSION.md').read_text();line=next(s for s in text.splitlines() if s.startswith('cd /d '));assert line.count(chr(92))==4,line
print('V63 final report, CMD syntax and preserved/new seals verified',flush=True)
