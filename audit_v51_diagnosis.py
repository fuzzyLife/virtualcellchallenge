"""Audit V51 diagnostic scope, reproducibility, summary arithmetic and frozen proofs."""
from pathlib import Path
import hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parent
A=ROOT/'analysis/v51_plateau'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 d=json.loads((A/'diagnosis.json').read_text());s=json.loads((A/'selection_sensitivity.json').read_text());t=json.loads((A/'rank_arithmetic_checks.json').read_text())
 assert len(d['cases'])==88
 assert len({(x['family'],x['case_index']) for x in d['cases']})==44
 assert len({(x['family'],x['target']) for x in d['cases']})==30
 assert len(d['cached_real_control_checks'])==5 and t['passed'] and t['fixtures_matched_scipy']==34
 assert s['leave_case_out'][0]['target']=='MTHFSD' and s['leave_case_out'][0]['leave_case_out_reach_delta']<0
 for summary in d['context_summary']:
  rows=[x for x in d['cases'] if x['family']==summary['family'] and x['arm']==summary['arm']]
  assert summary['cases']==len(rows)
  assert summary['median_true_DE']==np.median([x['true_DE'] for x in rows])
  assert summary['median_predicted_DE']==np.median([x['predicted_DE'] for x in rows])
 for panel in d['panel_sensitivity']:
  rows=[x for x in d['cases'] if x['family']==panel['family'] and x['arm']=='v24']
  assert panel['unique_targets']==len({x['target'] for x in rows})
 for item,script in [('diagnosis.json','diagnose_v51_plateau.py'),('selection_sensitivity.json','diagnose_v51_selection.py'),('rank_arithmetic_checks.json','test_v51_rank_claims.py')]:
  obj=json.loads((A/item).read_text());assert obj['script_sha256']==sha(ROOT/script)
 # Check every light-weight declared input independently; main diagnosis already
 # checked every population array against its original frozen array receipt.
 checked=0
 for obj in [d,s]:
  for name,h in obj['inputs_sha256'].items():
   p=ROOT/name
   if p.suffix=='.json':assert sha(p)==h,name;checked+=1
 pub=json.loads((ROOT/'experiments/v50_response_prior/analysis/published_result.json').read_text())
 for name,h in pub['frozen_pre_upload_receipt_sha256'].items():
  assert sha(ROOT/'experiments/v50_response_prior/analysis'/name)==h
 artifacts=[A/'diagnosis.json',A/'selection_sensitivity.json',A/'rank_arithmetic_checks.json',A/'case_diagnostics.csv',ROOT/'V51_PLATEAU_DIAGNOSIS.md']
 result=dict(passed=True,public_case_populations=44,distinct_context_target_pairs=30,arm_case_records=88,
  real_NTC_checks=5,rank_fixtures=34,metadata_hashes_rechecked=checked,frozen_preupload_receipts_unchanged=True,
  scope='Descriptive post-publication audit; no new model, blind test, oracle fit, release or leaderboard improvement claimed',
  script_sha256=sha(Path(__file__)),artifacts_sha256={str(p.relative_to(ROOT)):sha(p) for p in artifacts})
 (A/'audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
