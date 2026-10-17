"""Update README from actual V50 validation and delivery receipts."""
import json,re
import build_submission as b
def read(path):return json.loads(path.read_text()) if path.exists() else None
def main():
 A=b.ROOT/'experiments/v50_response_prior/analysis';decision=read(A/'release_decision.json');ready=read(A/'submission_ready.json');audit=read(A/'final_delivery_audit.json')
 status='The selected generic response-prior candidate passes development guards; additional-count confirmation is running.'
 if decision:status='Additional-count confirmation '+('passes' if decision['all_public_guards_passed'] else 'fails')+' the fixed public guards.'
 if ready:status+=' Official upload package is generated'+(' and independently verified.' if audit and audit['passed'] else '; final delivery audit pending.')
 section='\n'.join(['<!-- v38-v50:begin -->','','## V38–V50: independent evidence, population models and response priors (2026-10-17)','','Tested response activity/magnitude heads, KOLF/LINCS experimental evidence, contrastive identity objectives, population dispersion and ensembles, exact-margin projection, and bounded matrix rounding. V38–V49 do not establish a new release improvement; numerical tests alone did not translate to real-cell gains. V50 separates generic amplitude calibration from failed target-specific transfer. Both an averaged output-gene propensity prior and a fixed randomized checkpoint pass development guards; the frozen rule selects the randomized prior. Its discovery from a negative control is explicitly disclosed.','',status,'','[V45–V50 results, limitations and reproduction](V45_V50_RESULTS.md), [V38–V44 research ledger](V38_V44_RESULTS.md), [earlier methods](V38_V40_RESEARCH.md). V24 remains best published at0.14299785521745675. New private scores are unknown; D3 stays unopened. Scripts/receipts are versioned, large regenerable files are excluded. No automatic upload or Git push.','','<!-- v38-v50:end -->'])
 if ready and audit and audit['passed']:section=section.replace(status,status+' Upload file: [prediction.vcc](experiments/v50_response_prior/prediction.vcc).')
 path=b.ROOT/'README.md';text=path.read_text(encoding='utf-8-sig');pattern=r'<!-- v38-v(?:45|50):begin -->.*?<!-- v38-v(?:45|50):end -->'
 assert re.search(pattern,text,re.S);text=re.sub(pattern,lambda _:section,text,flags=re.S);path.write_text(text.rstrip()+'\n',encoding='utf-8')
if __name__=='__main__':main()

