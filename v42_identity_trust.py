"""V42: matched-objective contrast and source-validation abstention.

Subtract the calibration from the contrastive-label-shuffled control (whose
MSE auxiliary was correctly paired) from the genuine-identity calibration.
This is a matched-model contrast, not an exact gradient or causal decomposition.
"""
import argparse,json
import numpy as np
from threadpoolctl import threadpool_limits
import v40_identity_calibration as old
ROOT=old.ROOT;OUT=ROOT/'experiments/v42_identity_trust';A=OUT/'analysis'
OUTERS=old.OUTERS;ARMS=('contrast05','contrast10');data=old.data;b=old.b;emit=old.emit
def read(p):return old.read(p)
def write(p,v):old.write(p,v)
def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert b.sha256(ROOT/n)==h,n
 for n,h in p['inputs'].items():assert b.sha256(ROOT/n)==h,n
 return p
def chosen_validation(rec):
 return min(rec['history'],key=lambda h:h['objective'])
def calibration(outer):
 ident=old.load(outer,False);shuffled=old.load(outer,True);ri=old.read(old.A/'models'/(outer+'.json'));rs=old.read(old.A/'models'/(outer+'_shuffled.json'));vi=chosen_validation(ri);vs=chosen_validation(rs)
 use=vi['cross_entropy']<vs['cross_entropy'] and vi['cross_entropy']<ri['before']['cross_entropy']
 return ident-shuffled,{'outer':outer,'use_calibrator':bool(use),'identity_source_validation_CE':vi['cross_entropy'],'shuffled_source_validation_CE':vs['cross_entropy'],'baseline_source_validation_CE':ri['before']['cross_entropy'],'identity_weights_sha256':ri['sha256'],'shuffled_weights_sha256':rs['sha256'],'contrast_is_exact_gradient_decomposition':False}
def predict(delta,base,use=True,alpha=.10):
 return old.predict(delta,base,alpha) if use else base.copy()
def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v42_identity_trust.py','run_v42_identity_trust.py','v40_identity_calibration.py','run_v38_standard_counts.py','v38_response_hurdle.py','v16_methods.py','v12_cpu_metrics.py','v12_methods.py','v8_methods.py']
 inputs=['experiments/v40_identity_calibration/analysis/plan.json','experiments/v40_identity_calibration/analysis/selection.json','experiments/v41_identity_contrast/analysis/selection.json','experiments/v34_residual_flow/analysis/plan.json']
 for o in OUTERS:
  for stem in [o,o+'_shuffled']:
   for ext in ['npz','json']:inputs.append(str((old.A/'models'/(stem+'.'+ext)).relative_to(ROOT)))
 write(A/'plan.json',{'dependencies':{s:b.sha256(ROOT/s) for s in deps},'inputs':{s:b.sha256(ROOT/s) for s in inputs},'hypothesis':'After V41 coverage gains and amplitude/discrimination failures, test fixed5%/10% trust bounds; both evaluated with exact counts before selection. Retain only calibration specific to correctly labeled contrastive intervention matching; matched control retains same correctly paired MSE auxiliary. Source-validation CE must beat both baseline and matched control; otherwise exact V24 abstention. Difference between independently selected source-CV checkpoints is not an exact gradient decomposition.','coefficients':[.05,.10],'gain_bounds':{'contrast05':[.95,1.05],'contrast10':[.90,1.10]},'outer_exclusions':'Reuse exact V40 whole-study-excluded models. Source target-hash validation only determines abstention; no held-query outcomes used at inference.','primary':'All44 spent public means plus ideal pooled PDS, then original20 standard-generator cases with all six raw metrics. Frozen V38 count guards plus MSE/MAE<=1.01; at least two raw metrics improve; source objective beats shuffled when used. Null unchanged by zero identity. If passing, lock before additional24 count cases and secondseed380000; same count guards and per-study PDS/FID/REACH within0.03. No reselection after confirmation.','known_outcomes':'V38/V39/V40/V41 complete development outcomes and source validation controls known. Additional24 count cases and D3 reserve still unopened. Explicit outcome-informed next hypothesis on reused public data.','private_score_known':False})
 print('V42 confidence contrast frozen',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze']);a=p.parse_args();freeze()

