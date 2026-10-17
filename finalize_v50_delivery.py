"""Assemble a release manifest only after every independent delivery audit passes."""
from pathlib import Path
import json
import build_submission as b
ROOT=b.ROOT;A=ROOT/'experiments/v50_response_prior/analysis';OUT=A.parent
def read(n):return json.loads((A/n).read_text())
def main():
 ready=read('submission_ready.json');model=read('reference_model.json');decision=read('release_decision.json');validation=read('prediction_validation.json')
 assert ready['official_prep_passed'] and not ready['uploaded'] and decision['release'] and decision['all_public_guards_passed']
 assert (OUT/'prediction.vcc').stat().st_size==ready['bytes'];assert validation['n_cells']==360000 and validation['n_genes']==18533 and validation['nnz']<=4750000000
 audits=['preflight.json','independent_audit.json','raw_score_audit.json','production_inference_audit.json','source_row_audit.json','payload_roundtrip_audit.json','nullable_payload_contract_test.json','final_delivery_audit.json']
 for n in audits:assert read(n)['passed'],n
 assert read('independent_audit.json')['DE_set_comparisons']==171 and read('independent_audit.json')['profile_panels']==35
 assert read('source_row_audit.json')['rows']==360000 and read('production_inference_audit.json')['coordinates_replayed']==16679700
 basic=read('final_delivery_audit.json');assert basic['rows_independently_checked']==360000 and basic['groups']==900
 assert basic['submission_ready_sha256']==b.sha256(A/'submission_ready.json')
 assert basic['provenance_csv_sha256']==read('source_row_audit.json')['provenance_csv_sha256']
 assert read('preflight.json')['model_sha256']==ready['model_sha256']==b.sha256(A/'reference_model.json')
 assert read('production_inference_audit.json')['model_arrays_sha256']==model['arrays_sha256']
 assert ready['payload_auditor_sha256']==read('payload_roundtrip_audit.json')['auditor_sha256']==b.sha256(ROOT/'verify_v50_payload.py')
 payload=read('payload_roundtrip_audit.json')['full_CSR_comparison'];assert payload['data']==payload['indices']==validation['nnz'] and payload['indptr']==360001
 files=audits+['submission_ready.json','reference_model.json','release_decision.json','confirmation_stability.json','prediction_validation.json','runtime_receipt.json','official_runtime_receipt.json']
 b.write_json(A/'final_release_manifest.json',{'status':'Ready for manual validation-leaderboard upload; statistical and leaderboard superiority unproven','passed':True,'file':str(OUT/'prediction.vcc'),'bytes':ready['bytes'],'sha256':ready['sha256'],'selected_arm':model['selected_arm'],'prediction_shape':[360000,18533],'groups':900,'cells_per_group':400,'nnz':validation['nnz'],'aggregate_public_guards_passed':True,'seed_sensitive':True,'uploaded':False,'leaderboard_score':None,'best_published_reference':{'version':24,'overall':.14299785521745675},'proof_receipts_sha256':{n:b.sha256(A/n) for n in files},'finalizer_sha256':b.sha256(Path(__file__))})
 print('V50 final release manifest PASS',ready['bytes'],'bytes',flush=True)
if __name__=='__main__':main()


