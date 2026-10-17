"""Verify challenge controls are byte-identical to V24's production inputs."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v55_mean50_diagnostic/analysis'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
 return h.hexdigest()
def main():
 parent=ROOT/'experiments/v24_cd4_response/analysis/reference_model.json';d=json.loads(parent.read_text());rows=[]
 for r in d['context_original_inputs']:
  file=ROOT/r['file'];actual=sha(file);assert actual==r['sha256'],r['context'];rows.append({'context':r['context'],'file':r['file'],'bytes':file.stat().st_size,'sha256':actual});print('Original control hash matched',r['context'],flush=True)
 (A/'challenge_controls_audit.json').write_text(json.dumps({'passed':True,'records':rows,'parent_reference_model_sha256':sha(parent),'auditor_sha256':sha(Path(__file__))},indent=2)+'\n')
if __name__=='__main__':main()
