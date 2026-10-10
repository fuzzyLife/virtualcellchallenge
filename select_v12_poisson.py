"""Freeze a development winner before opening confirmation; never select from confirmation."""
import json
from pathlib import Path
from datetime import datetime,timezone
import pandas as pd
import build_submission as build
ROOT=build.ROOT;A=ROOT/'experiments/v12_distribution_calibrated/analysis/poisson_stage'
def main():
 destination=A/'count_selection.json';assert not destination.exists(),'Preserve an existing selection'
 guard=json.loads((A/'development_assessment.json').read_text());anchors=json.loads((A/'source_anchor_comparison.json').read_text());checks=json.loads((A/'source_checks.json').read_text())
 audit=json.loads((A/'development_evidence_audit.json').read_text());split=json.loads((A/'split_identity_audit.json').read_text())
 assert audit['passed'] and split['passed'] and split['plan_sha256']==build.sha256(A/'validation_plan.json')
 assert checks['passed'] and anchors['completed'];assert guard['scores_sha256']==build.sha256(A/'development_scores.csv')==anchors['scores_sha256']
 scores=pd.read_csv(A/'development_source_anchored_summary.csv').set_index('candidate');passing=[n for n in guard['passing_candidates'] if scores.loc[n,'overall']>scores.loc['v11','overall']]
 if not passing:
  build.write_json(A/'development_no_selection.json',{'created_utc':datetime.now(timezone.utc).isoformat(),'passing_raw_guards':guard['passing_candidates'],'passing_raw_and_anchored_improvement':[],'development_scores_sha256':build.sha256(A/'development_scores.csv'),'decision':'Retain V11; confirmation remains unopened.'});print('No candidate passes fixed raw guards and source-anchored improvement. V11 retained.');return
 name=max(passing,key=lambda n:float(scores.loc[n,'overall']))
 build.write_json(destination,{'created_utc':datetime.now(timezone.utc).isoformat(),'candidate':name,'confirmation_read':False,'development_scores_sha256':build.sha256(A/'development_scores.csv'),'plan_sha256':build.sha256(A/'validation_plan.json'),'methods_sha256':build.sha256(ROOT/'v11_methods.py'),'source_sha256':build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),'assessment_sha256':build.sha256(A/'development_assessment.json'),'anchors_sha256':build.sha256(A/'development_source_anchored_summary.csv'),'source_checks_sha256':build.sha256(A/'source_checks.json'),'script_sha256':build.sha256(Path(__file__)),'selection':'Fixed raw guards in both reference modes, then largest first-seed official public-source anchored mean among improvements. Private anchors unavailable; not a leaderboard forecast.'})
 print('Frozen development winner:',name,flush=True)
if __name__=='__main__':main()
