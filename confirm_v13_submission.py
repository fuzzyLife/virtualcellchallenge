"""Assess a frozen confirmation once; do not retune in response to its outcome."""
import json
from v13_methods import ROOT,A
from assess_v13_submission import assess
import assess_v13_submission as assessment_module
import build_submission as build

def main():
 destination=A/'confirmation_result.json';assert not destination.exists()
 choice=json.loads((A/'count_selection.json').read_text(encoding='utf-8'));folder=A/choice['validation_folder']
 assessment_module.A=folder;result=assess('confirmation')
 passed=choice['candidate'] in result['passing_candidates'] and choice['raw_guards_passed']
 build.write_json(destination,{'candidate':choice['candidate'],'confirmed_improvement':bool(passed),
  'development_raw_guards_passed':choice['raw_guards_passed'],'confirmation_raw_guards_passed':choice['candidate'] in result['passing_candidates'],
  'selection_sha256':build.sha256(A/'count_selection.json'),'scores_sha256':build.sha256(folder/'confirmation_scores.csv'),
  'summary_sha256':build.sha256(folder/'confirmation_summary.csv'),'assessment_sha256':build.sha256(folder/'confirmation_assessment.json'),
  'scope':'Public count-confirmation used once after development selection; source-summary biology has been viewed previously. Failure authorizes only an exploratory artifact, not a better-score claim.'})
 print('V13 confirmation recorded:',passed,flush=True)
if __name__=='__main__':main()
