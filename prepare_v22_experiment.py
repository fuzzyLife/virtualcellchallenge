"""Freeze V22 metadata-only development/oracle identities and release rules."""
import json,sys
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from v22_methods import ROOT,A,NAMES,write
import build_submission as b

def main():
 A.mkdir(parents=True,exist_ok=True);dest=A/'plan.json'
 if dest.exists():raise FileExistsError('Preserve frozen V22 plan')
 old=json.loads((ROOT/'experiments/v16_control_preserving/analysis/plan.json').read_text());contexts=[]
 for fi,parent in enumerate(old['contexts']):
  c={'family':parent['family'],'left':parent['left'],'right':parent['right'],'cells':[v for v in parent['cells'] if v['split']=='development'],'oracle':[]}
  names=sorted({v['target'] for v in c['cells']})[:3]
  for target in names:
   case=next(v for v in c['cells'] if v['target']==target);ids=np.random.default_rng(22001+fi*100+names.index(target)).permutation(case['actual']);template=np.random.default_rng(22002+fi).choice(c['left'],200,replace=False)
   c['oracle'].append({'target':target,'seed':22001+fi*100+names.index(target),'calibration':ids[:200].tolist(),'actual':ids[200:400].tolist(),'template':template.tolist()})
  contexts.append(c)
 deps=['v22_methods.py','prepare_v22_experiment.py','run_v22_experiment.py','assess_v22_experiment.py','collect_v22_cd4.py','extract_v22_cd4.py','test_v22_state_response.py','v16_methods.py','v8_methods.py','v6_methods.py','v12_cpu_metrics.py','v12_methods.py','fit_v7_response.py','audit_transfer_validation.py','build_submission.py','evaluate_atlas_transfer.py']
 plan={'version':22,'contexts':contexts,'names':NAMES,'eligible_finalists':['conditional'],'raw_guards':old['raw_guards'],'dependencies':{n:b.sha256(ROOT/n) for n in deps},'source_cache_sha256':old['source_cache_sha256'],'parent_plan_sha256':b.sha256(ROOT/'experiments/v16_control_preserving/analysis/plan.json'),
 'fit':'Training-only 768 stable non-panel state genes, four control-SVD factors; 2048 output genes. Matched source-batch controls; split-gene agreement regularizes continuous strengths. Separate shared vs target-specific fingerprints; ridge-160 state interactions, ridge-80 strength slopes; n/(n+200) and bounded .25 conditional LFC residual. Whole held family excluded from all fitting. Control reconstruction does not identify individual counterfactual states.',
 'oracle':'Nine old public target contrasts, disjoint 200 calibration/200 truth, five measured/predicted mean x compact/empirical-treated shape arms including actual calibration replica. Treated calibration depths shared across diagnostic arms. No 400+400 power claimed; no synthetic treatment duplicates.',
 'development':'14 previously used target/family pairs, two seeds, 400 true and predicted cells, seven arms. Count-baseline/fitting/source provenance audited. Not new independent biology.',
 'reserve':'New raw-cell CD4 D1 Rest panel selected alphabetically by metadata only, first eight competition targets with >=400 cells and >=40 matched controls per supported source batch, at most400 treatment cells per target. Disjoint control left/right halves within batch, >=16 each. Exclude all CD4 summary and raw data from deployable fitting. Same Marson study may supply prior summaries; not a new independent study. Reserve opened once only if development passes. If fewer than4 eligible targets, freeze explicit reduced-power amendment before scores or retain V16.',
 'release':'Existing available-metric raw guards + >=3 improved metrics including DE; centered response-projection Wasserstein improves; conditional beats scalar, independent and matched shuffled-state on >=2 of fidelity/reach/Jaccard. Null zero-effect identity exact V16. Only one finalist, same fixed comparisons/guards on raw CD4 reserve. Ineligible metrics remain missing; no official anchor score estimated. Failed candidates retain V16 and no improvement claim.',
 'fallback':'If the novel predictor is rejected, preserve champion V16. A requested upload-ready fallback may only be an explicitly labelled byte-identical V16 copy, never a new-model or improvement claim.'}
 write(dest,plan);print('Frozen V22 plan',b.sha256(dest),flush=True)
if __name__=='__main__':main()
