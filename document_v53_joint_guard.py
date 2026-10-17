"""Report V52 and feature-guarded V53 outcomes from saved evaluation receipts."""
import re
import v52_joint_response as m
from run_v53_joint_guard import A,BASE

def main():
 paths={'V52':BASE/'evaluation.json','V53 guarded':A/'evaluation.json'};summaries={}
 lines=['# Joint mean/rank predictor: V52 and corrected V53','',
  'Implemented and trained a shared 128/64 neural network for mean log2FC and rank effect relative to controls, with three whole-study outer exclusions, single-head ablations, linear heads and a trained shuffled-target control. All-source checkpoints use the same fixed fitting schedule. [V52 protocol](V52_JOINT_STRATEGY.md), [feature-contract correction](V53_JOINT_GUARD_STRATEGY.md), [original V52 results](V52_JOINT_RESULTS.md).','',
  'V53 reuses exact V52 labels and changes only how fitting-constant features are handled. It was frozen after partial predecessor evidence; it is not a blind new hypothesis. Both versions use repeatedly examined public biology; D3 remains unopened.','',
  '## Equal-study held-out errors','','| Version | Predictor | Mean MSE | Rank MSE |','|---|---|---:|---:|']
 for version,p in paths.items():
  if not p.exists():continue
  rec=m.read(p);summaries[version]=rec
  for arm,item in rec['summary'].items():
   x=item['equal_study'];lines.append(f"| {version} | {arm} | {x['mean_mse']:.7f} | {x['rank_mse']:.7f} |")
 if 'V53 guarded' in summaries:
  r=summaries['V53 guarded'];j=r['summary']['joint']['equal_study'];b=r['summary']['v24_pooled_rank']['equal_study'];z=r['summary']['zero_rank']['equal_study']
  lines+=['',f"Corrected joint mean MSE versus V24: **{100*(j['mean_mse']/b['mean_mse']-1):+.2f}%**. Rank MSE versus pooled rank: **{100*(j['rank_mse']/b['rank_mse']-1):+.2f}%**; versus zero rank: **{100*(j['rank_mse']/z['rank_mse']-1):+.2f}%**. Negative error changes are favorable.",'',
   'The zero-rank baseline retains V24 means and predicts theta=0; it is not a zero-expression submission. Mean/rank errors are raw response-label diagnostics, not official VCC scores. Directional and significance metrics require count construction and an independent reference.','',
   '## Prespecified pilot checks','']
  for k,v in r['checks'].items():lines.append(f'- {k}: **{v}**')
  lines+=['',f"Pilot success: **{r['pilot_success']}**.",'',
   ('The corrected joint model meets these response-level checks. A rank-aware count decoder and full population validation remain required before release.' if r['pilot_success'] else 'The corrected joint model does not meet the fixed response-level checks. Keep V24 as the preferred submission; the trained checkpoint is an exploratory artifact, not an established improvement.'),'',
   '## Study-specific errors','','| Held-out study | Joint mean / V24 MSE | Joint rank / pooled MSE | Joint rank / zero MSE |','|---|---:|---:|---:|']
  for group in m.OUTERS:
   ss=r['summary'];jj=ss['joint']['studies'][group];bb=ss['v24_pooled_rank']['studies'][group];zz=ss['zero_rank']['studies'][group]
   lines.append(f"| {group} | {jj['mean_mse']/bb['mean_mse']:.4f} | {jj['rank_mse']/bb['rank_mse']:.4f} | {jj['rank_mse']/zz['rank_mse']:.4f} |")
 lines+=['','## Label quality and limitations','',
  'There are 366 supported source-target groups: H1 12, HCT116 139, HEK293T 147, CD4 Rest D2 34 and CD4 Stim48 D2 34. All measured output genes are considered; fitting uses at most 768 eligible pairs per group and held-study evaluation uses all eligible measured outputs. Sample replicates are not new biological contexts. Hash0 targets are excluded from fitting labels, but allowed cross-study experimental responses to them remain retrieval features.','',
  'The paired-view rank correlations are approximately 0.691 in H1, 0.022 in HCT116, 0.249 in HEK293T and 0.62 in CD4. This weak Orion reproducibility makes noise shrinkage an especially important negative control. The complete source/view diagnostic is in V52 analysis/label_noise.json.','',
  'The feature guard prevents unsupported train-constant donor channels from activating untrained weights. It does not establish context-specific biological transfer or eliminate all distribution shift. Source mean labels and theta are lane-matched; query features use pooled control profiles. The model does not consume a fitted full control CDF/dispersion representation. A useful future extension would need source-validated control-distribution features and evidence that it beats zero/shuffled rank predictions, before a count decoder is trusted.','',
  '## Deliverables','',
  '- Corrected all-source checkpoint: experiments/v53_joint_guard/analysis/models/all/joint.pt, with SHA256 and fitting history in joint.json.',
  '- Corrected full-axis held-study mean/rank arrays: experiments/v53_joint_guard/analysis/predictions/{source}.npz, with source, masks and hashes in the JSON sidecars.',
  '- Original V52 models and arrays remain under experiments/v52_joint_mean_rank/analysis.',
  '- Root training, inference, tests and auditing scripts are versioned. Binary labels/models/arrays remain local and regenerable.','',
  '~~~python','from run_v53_joint_guard import install','install()','from predict_v52_joint_response import predict','result = predict(query_cpm, outer="all")','mean_effects = result["mean"]','rank_effects = result["rank"]','eligible = result["valid"]','~~~','',
  'query_cpm must contain 18,533 per-cell-mean CPM values in the competition gene order. The function returns all 300 targets by default; genes below the query control expression gate keep the V24 mean effect, with rank eligibility marked absent; the mask is not a guarantee of direct donor evidence. Own-target coordinates are excluded. The output is not an upload-ready count matrix.','',
  '## Audit status','']
 for name in ['label_audit.json','boundary_audit.json','feature_contract_audit.json','final_audit.json','delivery_manifest.json']:
  p=A/name;lines.append('- '+name+': '+('**passed**' if p.exists() and m.read(p)['passed'] else 'pending'))
 lines+=['','The [delivery manifest](experiments/v53_joint_guard/analysis/delivery_manifest.json) binds the all-source checkpoint and API code to completed audits and a CPU reload/forward smoke test. No new .vcc file was created or uploaded. No Git push was performed.','']
 (m.ROOT/'V53_JOINT_GUARD_RESULTS.md').write_text('\n'.join(lines),encoding='utf-8')
 status='Corrected retraining is in progress.'
 if 'V53 guarded' in summaries:
  r=summaries['V53 guarded'];status='Training and whole-study evaluation completed. Fixed response-level pilot checks '+('pass' if r['pilot_success'] else '**fail**')+'.'
 section='<!-- v53-joint:begin -->\n\n## V53 corrected joint mean/rank predictor\n\n'+status+' Preserves V52 and reuses exact paired-view labels, disabling train-constant donor channels at fitting and inference. Architecture, settings, folds and controls are unchanged. [Results and model API](V53_JOINT_GUARD_RESULTS.md), [correction protocol](V53_JOINT_GUARD_STRATEGY.md). No count submission or private-score improvement is claimed; D3 remains unopened.\n\n<!-- v53-joint:end -->'
 p=m.ROOT/'README.md';text=p.read_text(encoding='utf-8-sig');pattern=r'<!-- v53-joint:begin -->.*?<!-- v53-joint:end -->'
 text=re.sub(pattern,lambda _:section,text,flags=re.S) if re.search(pattern,text,re.S) else text.replace('# Virtual Cell Challenge 2026','# Virtual Cell Challenge 2026\n\n'+section,1);p.write_text(text,encoding='utf-8')
 print('Joint predictor report updated',flush=True)
if __name__=='__main__':main()
