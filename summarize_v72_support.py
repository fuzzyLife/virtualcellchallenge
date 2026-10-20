"""Compare nominal and actual response-evidence coverage before/after batch pooling."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
A=ROOT/'experiments/v72_batch_offset/analysis'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 rows=[];inputs={}
 for p in sorted((A/'sources').glob('*.json')):
  old=ROOT/'experiments/v70_count_teacher/analysis/sources'/p.name
  l=json.loads(old.read_text());r=json.loads(p.read_text());left={x['target']:x for x in l['targets']};right={x['target']:x for x in r['targets']}
  assert set(left)==set(right),p.name
  assert all(right[t]['cells']==left[t]['cells'] for t in left)
  a=sum(x['supported_genes']>0 for x in left.values());b=sum(x['supported_genes']>0 for x in right.values())
  assert b>=a
  rows.append(dict(source=p.stem,nominal_targets=len(left),V70_supported_targets=a,V72_supported_targets=b,targets_with_new_evidence=b-a))
  inputs[old.relative_to(ROOT).as_posix()]=sha(old);inputs[p.relative_to(ROOT).as_posix()]=sha(p)
 assert len(rows)==6
 out=dict(rows=rows,inputs_sha256=inputs,script_sha256=sha(Path(__file__)),same_treated_cell_counts=True,explanation='Pooling matched batch offsets before effect estimation removes the >=10-cells-per-batch censor. Eligibility still requires >=40 total matched treated cells.')
 (A/'coverage_comparison.json').write_bytes((json.dumps(out,indent=2)+'\n').encode())
 print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
