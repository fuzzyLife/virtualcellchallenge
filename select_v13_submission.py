"""Freeze one model before opening count confirmation; preserve all development failures."""
import json
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd
from v13_methods import ROOT,A
import build_submission as build

def main():
 assert not (A/'count_selection.json').exists()
 replay=json.loads((A/'baseline_replay_audit.json').read_text(encoding='utf-8'));assert replay['passed']
 records=[]
 for folder,module,checksfile in [(A,'v13_methods','response_checks.json'),(A/'robust_stage','v13_robust_methods','robust_source_checks.json'),(A/'quality_stage','v13_quality_methods','quality_source_checks.json')]:
  assessment=json.loads((folder/'development_assessment.json').read_text(encoding='utf-8'))
  anchor=json.loads((folder/'source_anchor_comparison.json').read_text(encoding='utf-8'));assert anchor['completed']
  checks=json.loads((A/checksfile).read_text(encoding='utf-8'));assert checks['passed']
  summary=pd.read_csv(folder/'development_summary.csv');scaled=pd.read_csv(folder/'development_source_anchored_summary.csv').set_index('candidate')
  base=summary[(summary.candidate=='v11')&summary.reference_mode.eq('full_independent')].iloc[0]
  for _,row in summary[summary.candidate.ne('v11')&summary.reference_mode.eq('full_independent')].iterrows():
   name=row.candidate
   records.append({'candidate':name,'validation_folder':str(folder.relative_to(A)) if folder!=A else '',
    'methods_module':module,'source_checks_file':checksfile,'raw_guards_passed':name in assessment['passing_candidates'],
    'expression_reduction':float(base.common_expression_ratio-row.common_expression_ratio),
    'source_anchored_overall':float(scaled.loc[name,'overall']),'baseline_source_anchored_overall':float(scaled.loc['v11','overall'])})
 eligible=[r for r in records if r['raw_guards_passed']]
 if eligible:chosen=max(eligible,key=lambda r:r['expression_reduction']);policy='Prespecified raw guards, then greatest primary expression reduction.'
 else:chosen=max(records,key=lambda r:r['source_anchored_overall']-r['baseline_source_anchored_overall']);policy='No raw-guard winner: freeze the largest public anchored score change for an explicitly exploratory artifact. No confirmed-gain claim.'
 folder=A/chosen['validation_folder'];selection={'created_utc':datetime.now(timezone.utc).isoformat(),**chosen,'confirmation_read':False,'development_comparisons':records,
  'development_scores_sha256':build.sha256(folder/'development_scores.csv'),'plan_sha256':build.sha256(folder/'validation_plan.json'),
  'methods_sha256':build.sha256(ROOT/(chosen['methods_module']+'.py')),'selection_policy':policy,
  'baseline_replay_sha256':build.sha256(A/'baseline_replay_audit.json'),'script_sha256':build.sha256(Path(__file__))}
 build.write_json(A/'count_selection.json',selection)
 if folder!=A:build.write_json(folder/'count_selection.json',selection)
 print('Frozen V13 candidate:',chosen,'policy:',policy,flush=True)
if __name__=='__main__':main()
