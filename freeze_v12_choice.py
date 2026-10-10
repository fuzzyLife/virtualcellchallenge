"""Choose once across the completed bounded development stages, before any confirmation."""
import json,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd
import build_submission as build
ROOT=build.ROOT;A=ROOT/'experiments/v12_distribution_calibrated/analysis'
STAGES={'flex_stage':'select_v12_flex.py','harmonized_stage':'select_v12_harmonized.py','poisson_stage':'select_v12_poisson.py'}
def main():
 destination=A/'final_development_choice.json';assert not destination.exists()
 choices=[];evidence={}
 for folder,selector in STAGES.items():
  d=A/folder;assert not (d/'confirmation_scores.csv').exists() and not (d/'count_selection.json').exists()
  if (d/'futility_stop.json').exists():
   stopped=json.loads((d/'futility_stop.json').read_text());assert stopped['terminal_rejection'] and all(v['terminal_rejection'] for v in stopped['results'].values())
   assert stopped['plan_sha256']==build.sha256(d/'validation_plan.json') and stopped['scores_sha256']==build.sha256(d/'development_scores.csv')
   evidence[folder]={name:build.sha256(d/name) for name in ['validation_plan.json','futility_stop.json','source_checks.json','split_identity_audit.json']};continue
  guard=json.loads((d/'development_assessment.json').read_text());audit=json.loads((d/'development_evidence_audit.json').read_text());anchor=json.loads((d/'source_anchor_comparison.json').read_text());assert audit['passed'] and anchor['completed']
  assert guard['scores_sha256']==anchor['scores_sha256']==build.sha256(d/'development_scores.csv')
  scores=pd.read_csv(d/'development_source_anchored_summary.csv').set_index('candidate');base=float(scores.loc['v11','overall'])
  for candidate in guard['passing_candidates']:
   gain=float(scores.loc[candidate,'overall'])-base
   if gain>0:choices.append({'stage':folder,'candidate':candidate,'development_source_anchored_gain':gain})
  evidence[folder]={name:build.sha256(d/name) for name in ['validation_plan.json','development_assessment.json','development_evidence_audit.json','development_source_anchored_summary.csv','source_checks.json','split_identity_audit.json']}
 if choices:
  choice=max(choices,key=lambda c:c['development_source_anchored_gain']);subprocess.run([sys.executable,str(ROOT/STAGES[choice['stage']])],check=True)
  selected=A/choice['stage']/'count_selection.json';assert json.loads(selected.read_text())['candidate']==choice['candidate'];choice['selection_sha256']=build.sha256(selected)
 else:choice={'stage':None,'candidate':'v11','decision':'No new candidate passes fixed raw guards and public source-anchored improvement; retain V11, leave confirmation unopened.'}
 build.write_json(destination,{'created_utc':datetime.now(timezone.utc).isoformat(),'choice':choice,'passing_choices':choices,'evidence':evidence,'confirmation_read':False,'script_sha256':build.sha256(Path(__file__)),'scope':'One final development choice across three bounded source/emitter stages, all on reused public development panels. Confirmation may be used once; no retuning after failure. Source anchors are noisy and are not private leaderboard anchors.'})
 print(json.dumps(choice,indent=2),flush=True)
if __name__=='__main__':main()
