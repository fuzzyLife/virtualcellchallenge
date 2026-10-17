"""Render completed V52 training/evaluation status without inventing missing stages."""
import json,re
import v52_joint_response as m

def main():
 lines=['# V52 joint mean/rank training results','',
  'Shared neural predictor with separate mean log2FC and Mann-Whitney rank-effect outputs. [Frozen protocol and reproduction](V52_JOINT_STRATEGY.md). Public data are reused; D3 remains unopened. Mean/rank arrays require separate count-decoder validation before a competition release.','',
  '## Extracted paired-view labels','','| Source | Study | Supported targets | Eligible target/gene pairs |','|---|---|---:|---:|']
 for source in m.SOURCES:
  p=m.A/'labels'/(source+'.json')
  if p.exists():
   r=m.read(p);lines.append(f"| {source} | {m.GROUP[source]} | {len(r['groups'])} | {sum(x['valid_genes'] for x in r['groups']):,} |")
  else:lines.append(f'| {source} | {m.GROUP[source]} | pending | pending |')
 lines+=['','The two views are sampling replicates, not independent studies. Missing target/gene labels abstain. All measured genes were considered; each fitting target samples at most 768 eligible gene pairs, while outer evaluation uses every eligible measured gene.','',
  '## Fitted models','','| Excluded study | Model | Selected epoch | Training pairs | Source-validation pairs |','|---|---|---:|---:|---:|']
 for outer in [*m.OUTERS,'all']:
  for arm in (('joint',) if outer=='all' else m.ARMS):
   p=m.A/'models'/outer/(arm+'.json')
   if p.exists():
    r=m.read(p);lines.append(f"| {outer} | {arm} | {r['best_epoch']} | {r['trained_pairs']:,} | {r['validation_pairs']:,} |")
 p=m.A/'evaluation.json';summary=None
 if p.exists():
  evaluation=m.read(p);summary=evaluation['summary'];lines+=['','## Whole-study raw prediction errors','','| Predictor | Mean MSE | Rank MSE | Mean MAE | Rank MAE |','|---|---:|---:|---:|---:|']
  for arm,item in summary.items():
   v=item['equal_study'];lines.append(f"| {arm} | {v['mean_mse']:.6f} | {v['rank_mse']:.6f} | {v['mean_mae']:.6f} | {v['rank_mae']:.6f} |")
  lines+=['','Prespecified pilot checks:']
  for k,v in evaluation['checks'].items():lines.append(f'- {k}: **{v}**')
  lines+=['',f"Pilot success: **{evaluation['pilot_success']}**. This is not a private score or an upload decision.",'','## Per-study joint versus baseline','', '| Study | Joint mean MSE / V24 | Joint rank MSE / pooled | Joint rank MSE / zero |','|---|---:|---:|---:|']
  for group in m.OUTERS:
   j=summary['joint']['studies'][group];b=summary['v24_pooled_rank']['studies'][group];z=summary['zero_rank']['studies'][group]
   lines.append(f"| {group} | {j['mean_mse']/b['mean_mse']:.4f} | {j['rank_mse']/b['rank_mse']:.4f} | {j['rank_mse']/z['rank_mse']:.4f} |")
  joint=summary['joint']['equal_study'];base=summary['v24_pooled_rank']['equal_study']
  lines+=['',f"Joint changes versus the reference: mean MSE {100*(joint['mean_mse']/base['mean_mse']-1):+.2f}%, rank MSE {100*(joint['rank_mse']/base['rank_mse']-1):+.2f}%. Negative error changes are favorable. View disagreement, few independent studies, assay differences and the bounded residual architecture limit the interpretation. The same raw prediction loss does not ensure better significance, FID, REACH or a private rank."]
 else:lines+=['','Whole-study evaluation is pending.']
 lines+=['','## Integrity and delivery','']
 for n in ['label_audit.json','boundary_audit.json','final_audit.json']:
  p=m.A/n;lines.append(f"- {n}: "+('**passed**' if p.exists() and m.read(p)['passed'] else 'pending'))
 lines+=['','Model files live under experiments/v52_joint_mean_rank/analysis/models. The all-source joint.pt checkpoint is an exploratory inference artifact; it does not have an independent all-source test. Full-axis held-study joint mean/rank predictions are under analysis/predictions. Each model and prediction matrix has a SHA256 sidecar. Large binaries are ignored by Git.','',
  'Root API: predict_v52_joint_response.predict(query_cpm, outer="all") returns mean effects, rank effects and a supported-coordinate mask. It keeps V24 mean effects for output genes below the rank-label expression gate and excludes the own gene. A rank effect is not a p-value. Current source-control features are means/ranks and donor response statistics, not a fitted full query-control distribution; a later rank-aware count decoder would need that distribution and tie structure. Production inference has more available donor studies than some fitting folds; this requires further count and biological validation.','',
  'Audits: independent raw-cell pairwise rank replay; treatment/control view disjointness; excluded-response mutation; source-validation target-group boundaries; model and inference array hashes/bounds. These checks verify the implementation and exclusions, not hidden-context generalization. No new .vcc file was created or submitted.','',
  '~~~bat','SCIENCE_PYTHON audit_v52_joint_response.py labels','SCIENCE_PYTHON audit_v52_joint_response.py boundaries','SCIENCE_PYTHON audit_v52_joint_response.py artifacts','SCIENCE_PYTHON document_v52_joint_response.py','~~~','']
 (m.ROOT/'V52_JOINT_RESULTS.md').write_text('\n'.join(lines),encoding='utf-8')
 status='Extraction/training in progress.' if not summary else 'Training and whole-study raw prediction evaluation completed; see the measured results and pilot checks.'
 section='<!-- v52-joint:begin -->\n\n## V52 joint mean/rank predictor\n\n'+status+' Shared 128/64 neural network, single-head ablations, independent linear heads and a trained shuffled-target control. Whole H1/Orion/CD4 study exclusions; two disjoint raw population views supply full-axis mean/rank labels. [Protocol](V52_JOINT_STRATEGY.md), [results and reproducibility](V52_JOINT_RESULTS.md). This is model training and response-array evaluation, not a validated count submission. D3 remains unopened.\n\n<!-- v52-joint:end -->'
 path=m.ROOT/'README.md';text=path.read_text(encoding='utf-8-sig');pattern=r'<!-- v52-joint:begin -->.*?<!-- v52-joint:end -->'
 text=re.sub(pattern,lambda _:section,text,flags=re.S) if re.search(pattern,text,re.S) else text.replace('# Virtual Cell Challenge 2026','# Virtual Cell Challenge 2026\n\n'+section,1)
 path.write_text(text,encoding='utf-8')
if __name__=='__main__':main()
