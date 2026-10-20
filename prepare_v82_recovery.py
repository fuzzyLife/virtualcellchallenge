"""Prepare a byte-identical champion recovery only after all new pilots fail.

This is an upload artifact, not a newly learned predictor or new score claim.
Uses the official CLI environment and preserves every prior artifact.
"""
from pathlib import Path
import hashlib,json,shutil
from datetime import datetime,timezone
from vcc.vccfile import validate_vcc,read_vcc_meta
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v82_v66_recovery';A=OUT/'analysis'

def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(16*1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(p,v):p.write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')

def rejected():
 records=[]
 for v in range(70,82):
  folder=next((ROOT/'experiments').glob(f'v{v}_*'))/'analysis';selection=read(folder/'selection.json');proof=[folder/'selection.json']
  if v==75:
   p=folder/'broader_assessment.json';r=read(p);assert not r['passes'];proof.append(p)
  elif v==76:
   p=folder/'validation_complete.json';r=read(p);assert not r['primary_critical_gates_passed'] and not r['broad_passed'];proof.append(p)
  else:assert selection['selected'] is None,(v,'A new candidate qualified; evaluate it before any recovery release')
  if v>=78:
   p=folder/'development_0_independent_audit.json';assert read(p)['passed'];proof.append(p)
  assert not (folder.parent/'prediction.vcc').exists(),(v,'Existing new-model artifact needs review')
  records.append(dict(version=v,selected=selection.get('selected'),qualified_for_release=False,evidence={str(p.relative_to(ROOT)):digest(p) for p in proof}))
 return records

def main():
 evidence=rejected();source=ROOT/'experiments/v66_assay_strength/prediction.vcc';original=ROOT/'experiments/v66_assay_strength/analysis';audit=read(original/'independent_validation.json');receipt=read(original/'leaderboard_result.json');expected=audit['package_sha256'];assert audit['passed'] and expected==receipt['package_sha256']
 A.mkdir(parents=True,exist_ok=True);destination=OUT/'prediction.vcc';before=source.stat()
 if not destination.exists():
  pending=destination.with_suffix('.vcc.partial');assert not pending.exists();h=hashlib.sha256();n=0
  with source.open('rb') as src,pending.open('xb') as dst:
   for b in iter(lambda:src.read(16*1024*1024),b''):h.update(b);dst.write(b);n+=len(b)
  after=source.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
  assert h.hexdigest()==expected and n==receipt['upload_bytes'];assert digest(pending)==expected;pending.replace(destination)
 else:assert digest(destination)==expected
 validate_vcc(destination);meta=read_vcc_meta(destination);assert meta['nnz']==audit['packaged_nonzero_values_checked']
 record=dict(file=str(destination),sha256=expected,bytes=destination.stat().st_size,official_container_validated=True,byte_identical_to_V66=True,new_predictor=False,new_hidden_score=False,uploaded=False,created_utc=datetime.now(timezone.utc).isoformat(),parent_entry=receipt['entry_id'],parent_published_overall=receipt['overall'],parent_rank_is_historical=True,parent_full_payload_audit_sha256=digest(original/'independent_validation.json'),prior_public_guards_passed=False,new_methods=evidence,package_meta=meta,script_sha256=digest(Path(__file__)),interpretation='Recovery of best published project model after V69 regression. No new model met its frozen public guards. This copy does not establish improvement beyond V66 and resubmitting consumes an allowance.')
 write(A/'submission_ready.json',record);print(json.dumps({k:record[k] for k in ['file','sha256','bytes','byte_identical_to_V66','new_predictor','uploaded']},indent=2))
if __name__=='__main__':main()
