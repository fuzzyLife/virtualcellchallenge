"""Verify Git preserves the bytes of all newly staged scientific artifacts."""
import hashlib,json,subprocess
from pathlib import Path
import build_submission as b
ROOT=b.ROOT
def main():
 base=['git','-c','safe.directory='+str(ROOT).replace('\\','/')]
 def run(args,**kw):return subprocess.run(base+args,cwd=ROOT,check=True,stdout=subprocess.PIPE,**kw).stdout
 names=[n.decode() for n in run(['diff','--cached','--name-only','-z']).split(b'\0') if n]
 names=[n for n in names if n not in ['README.md','.gitignore','.gitattributes'] and n!='experiments/v45_bounded_dispersion/analysis/staged_byte_audit.json']
 payload=run(['cat-file','--batch'],input=b''.join((':'+n+'\n').encode() for n in names));offset=0;checks=[]
 for name in names:
  end=payload.index(b'\n',offset);head=payload[offset:end].decode();size=int(head.split()[-1]);data=payload[end+1:end+1+size];offset=end+size+2
  file=ROOT/name;assert data==file.read_bytes(),'Staged bytes differ: '+name
  checks.append({'file':name,'sha256':hashlib.sha256(data).hexdigest()})
 receipt={'passed':True,'files':len(checks),'checks':checks,'excluded':'README and Git configuration have conventional EOL normalization; audit receipt excludes itself. Scientific scripts and frozen metadata require exact working/index bytes.','auditor_sha256':b.sha256(Path(__file__))}
 b.write_json(ROOT/'experiments/v45_bounded_dispersion/analysis/staged_byte_audit.json',receipt);print('Staged scientific byte audit PASS',len(checks),flush=True)
if __name__=='__main__':main()

