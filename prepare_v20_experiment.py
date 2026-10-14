"""Freeze Stack comparison only after source preparation and author-output parity pass."""
import json
from pathlib import Path
import build_submission as b
from v20_methods import ROOT,A,R,FAMILIES,NAMES

def main():
 A.mkdir(parents=True,exist_ok=True)
 if (A/'plan.json').exists():raise FileExistsError('Preserve frozen V20 comparison')
 parent=ROOT/'experiments/v19_covariance_forces/analysis';old=json.loads((parent/'plan.json').read_text());assert not json.loads((parent/'selection.json').read_text())['new_candidate_passes'];parity=json.loads((A/'mean_adapter_parity.json').read_text());assert parity['passed'] and parity['adapter_sha256']==b.sha256(ROOT/'stack_v20_mean_adapter.py');prep=json.loads((A/'development_input_preparation.json').read_text());assert prep['complete']
 deps=['v20_methods.py','stack_v20_mean_adapter.py','prepare_v20_inputs.py','prepare_v20_experiment.py','infer_v20_stack.py','run_v20_experiment.py','assess_v20_experiment.py','audit_v20_experiment.py','test_v20_counts.py','v19_methods.py','v16_methods.py','fit_v7_response.py','v12_cpu_metrics.py','v8_methods.py','v6_methods.py','build_submission.py','evaluate_atlas_transfer.py'];files=[]
 for family in FAMILIES:
  for suffix in ['.json','_query.npz','_counts.npz']:
   path=A/'inputs'/(family+suffix);files.append({'path':str(path.relative_to(ROOT)).replace('\\','/'),'sha256':b.sha256(path)})
 for family,ti in prep['source_target_pairs']:
  import pandas as pd
  targets=pd.read_csv(ROOT/'pert_counts.csv').iloc[:,0].astype(str).tolist();target=targets[ti]
  for suffix in ['.json','_control.npz','_treated.npz']:
   path=A/'inputs/sources'/(family+'_'+target+suffix);files.append({'path':str(path.relative_to(ROOT)).replace('\\','/'),'sha256':b.sha256(path)})
 record={'version':20,'contexts':old['contexts'],'seeds':old['seeds'],'family_order':FAMILIES,'candidates':NAMES,'eligible':['stack_mean','stack_joint'],'raw_guards':old['raw_guards'],'dependencies':{n:b.sha256(ROOT/n) for n in deps},'input_files':files,'source_cache_sha256':old['source_cache_sha256'],'parent_plan_sha256':b.sha256(parent/'plan.json'),'previous_failed_selection_sha256':b.sha256(parent/'selection.json'),'assets_manifest_sha256':b.sha256(A/'assets.json'),'checkpoint_sha256':'f93cf6f42f36c8a85dc570d92e801c1fc3e1f45d55741e6bb250d178a6b6ad36','checkpoint_revision':'b09f085dac03d170b078a5c72f550ae93686e544','author_code_revision':'cacc2e4b09435c3e536d46237d10b50f222dd144','mean_adapter_parity_sha256':b.sha256(A/'mean_adapter_parity.json'),
 'decoder':'Original pretrained aligned model, strict full weights, FP32, batch4, one step with prompt .25/context .4; exactly original get_prediction mean/dispersion outputs (bitwise certified). Sampler/MDM iterations omitted, no approximation to weights. This is a bounded one-step expert, not reproduction of the authors five-step drug demonstration. Final counts use existing validated integer generation.',
 'source_prompts':'Whole query family excluded from prompts. Choose available >=40-source-cell family by direct pairwise valid Spearman of NTC mean CPM, all panel genes excluded, no +1/log transform. Query controls only; matched source NTC/trt batches row by row. Same family target-shuffle within [40,80),[80,200),[200,infinity) source-size bins; source sizes matched by prefix min(real,shuffled,400). Duplicated source NTC conditioning rows disclosed, not counted as independent treatments.',
 'mean':'Within common measured output vocabulary, calibrate log2(meanCPM generated treatment / generated NTC). Variance combines conditional-mean cell variance and analytic per-query NB mu+mu^2/theta, with independent-arm approximation. Weight=(lfc^2/(lfc^2+variance+.0025))*.04/(.04+variance)*n/(n+100). Update V16 by .25*weight*clip(StackLFC-V16,-.5,.5) plus original expression gate; <=.125 log2 correction. Missing outputs retain V16 multipliers. Uncalibrated generated treatment vs real NTC is a diagnostic, never eligible.',
 'joint':'Paired conditional mean log1p CP10k shifts, centered per gene, projected into existing training-only16 control programs. Map cached query states to each template by nearest program scores. Decode residual via program basis/query derivative, clip +/-1 log2 and weight.25; modify only the .10 empirical component of residual10. Exact mean/bulk count constraints and count seeds retained. Source n/(n+100) also weights residual. Zero residual routes bitwise to V16.',
 'reweight':'Same calibrated means, but tilt real query states along generated mean program direction; target-dependent bounded exponential weights and seeded empirical-component resampling, keeping original output libraries. This is a simpler occupancy control; not identified cell birth/death. Amplitude comparator matches Stack mean non-own L1. Shuffled prompts use the same support/bounds.',
 'shape':'Center each population in training16-program log CP10k space; mean Wasserstein along32 fixed unit Gaussian directions(seed20020), relative covariance Frobenius error and measured nonpanel detection MAE. Scoring uses actual aligned count library totals, while input generation retains full-source totals where needed.',
 'release':'Available raw/worst-family guards, >=3 improved metrics including DE, beat shuffled and amplitude controls on >=2 fidelity/reach/Jaccard. Joint additionally beats mean and empirical reweighting controls on >=2 DE measures and lower centered shape Wasserstein. Same once-only rule on unspent8 reserved pairs; missing HCT NMAE remains missing. No leaderboard guarantee.',
 'pretraining_overlap':'NOT independently audited. Stack pretrained on scBaseCount150M and aligned CXG45M/Parse10M. Prompts exclude query family; pretrained weights may have encountered public benchmark studies. Do not describe this as a clean unseen-pretraining-family test.',
 'verification':'Corrected V17 eight-pair algorithm reserve, not opened by rejected V17/V18/V19. Old accessible biology partly reused; not blind independent biological confirmation.'}
 b.write_json(A/'plan.json',record);print('V20 frozen',b.sha256(A/'plan.json'),flush=True)
if __name__=='__main__':main()
