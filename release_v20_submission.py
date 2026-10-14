"""Declare V20 upload-ready only after complete independent validation."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v20_stack_transfer/analysis';OUT=A.parent

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
 return h.hexdigest()
def main():
 ready=json.loads((A/'submission_ready.json').read_text());independent=json.loads((A/'independent_validation.json').read_text());inputs=json.loads((A/'production_input_independent_audit.json').read_text());confirmation=json.loads((A/'confirmation_result.json').read_text());assert confirmation['new_candidate_passes'] and confirmation['candidate']=='stack_mean';assert independent['passed'] and inputs['passed'];assert independent['verified_rows']==360000 and independent['verified_groups']==900 and independent['exact_replayed_groups']==12;assert independent['all_production_responses_independently_reconstructed'] and inputs['query_controls_verified']==1200;assert independent['package_sha256']==ready['sha256']==sha(OUT/'prediction.vcc');assert independent['model_sha256']==sha(A/'reference_model.json');assert inputs['production_plan_sha256']==sha(A/'production_plan.json')
 compatibility=json.loads((A/'production_runtime_compatibility.json').read_text());assert compatibility['wrapper_sha256']==sha(ROOT/'run_v20_submission_inference.py') and compatibility['frozen_inference_sha256']==sha(ROOT/'infer_v20_submission.py') and compatibility['pandas_future_infer_string'] is False
 ready.update(inference_runtime_compatibility_sha256=sha(A/'production_runtime_compatibility.json'),independent_full_payload_validation_passed=True,independent_production_input_validation_passed=True,upload_ready=True,independent_validation_sha256=sha(A/'independent_validation.json'),production_input_audit_sha256=sha(A/'production_input_independent_audit.json'),release_checker_sha256=sha(Path(__file__)),uploaded=False,leaderboard_score=None)
 (A/'submission_ready.json').write_text(json.dumps(ready,indent=2)+'\n',encoding='utf-8');print('V20 is independently validated and ready for manual upload:',OUT/'prediction.vcc',flush=True)
if __name__=='__main__':main()
