"""Assess frozen V9 once without confirmation-based reselection."""
import json
from v9_methods import ROOT,A
from select_v9_generator import assess
import build_submission as build

if __name__=='__main__':
    if (A/'confirmation_result.json').exists():raise FileExistsError('Already confirmed')
    selection=json.loads((A/'count_selection.json').read_text());winner=selection['candidate']
    assert selection['methods_sha256']==build.sha256(ROOT/'v9_methods.py')
    assert selection['evaluation_script_sha256']==build.sha256(ROOT/'evaluate_v9_submission.py')
    score,fidelity,guards=assess('confirmation');difference=float(score.loc[winner,'utility']-score.loc['legacy','utility'])
    build.write_json(A/'confirmation_result.json',{'candidate':winner,'guard_passed':guards[winner],
        'confirmed_improvement':bool(guards[winner] and difference>0),'utility_difference':difference,
        'supported_budget_fidelity':fidelity,'summary_sha256':build.sha256(A/'count_confirmation_summary.csv'),
        'frozen_calibration_selection_sha256':build.sha256(A/'count_selection.json'),
        'policy':'Frozen calibration winner is assessed once; no reselection or parameter adjustment from confirmation.'})
    print('V9 confirmation',winner,difference,guards[winner],flush=True)
