"""Apply an independent neutral-control gate to the frozen V8 calibration candidates."""
import json
import pandas as pd
import numpy as np
import build_submission as build
from v8_methods import ROOT,A


def supported_budget_fidelity(frame):
    supported=frame[frame.reference_calls>=10]
    result={}
    for candidate,rows in supported.groupby('candidate'):
        budget=np.maximum(rows.reference_calls,rows.predicted_calls)
        result[candidate]=float((rows.fidelity.fillna(0)*budget).sum()/budget.sum())
    return result


def main():
    if (A/'safe_count_selection.json').exists():raise FileExistsError('Safe calibration winner is already frozen')
    if (A/'confirmation_result.json').exists():raise ValueError('Safety selection must precede confirmation scoring')
    original=json.loads((A/'count_selection.json').read_text());neutral=json.loads((A/'decoder_neutral_manifest.json').read_text())
    assert original['calibration_summary_sha256']==build.sha256(A/'count_calibration_summary.csv')
    assert neutral['scores_sha256']==build.sha256(A/'decoder_neutral_scores.csv')
    assert neutral['methods_sha256']==build.sha256(ROOT/'v8_methods.py')
    score=pd.read_csv(A/'count_calibration_summary.csv').set_index('candidate')
    cell_scores=pd.read_csv(A/'count_calibration_scores.csv')
    fidelity=supported_budget_fidelity(cell_scores);base=score.loc['baseline_legacy']
    allowed=[]
    for name,setting in original['settings'].items():
        row=score.loc[name]
        stable_guard=fidelity.get(name,0)>=fidelity.get('baseline_legacy',0)-.02
        if (setting[1] in neutral['passed_generators'] and stable_guard and
            row.common_expression_ratio<=1.02*base.common_expression_ratio and row.pds>=base.pds-.05):allowed.append(name)
    if not allowed:raise ValueError('No calibrated candidate passes the independent neutral-control gate')
    winner=max(allowed,key=lambda n:score.loc[n,'utility'])
    record={**original,'candidate':winner,'response':original['settings'][winner][0],'generator':original['settings'][winner][1],
            'passing_safe_candidates':allowed,'original_count_winner':original['candidate'],
            'neutral_manifest_sha256':build.sha256(A/'decoder_neutral_manifest.json'),
            'original_count_selection_sha256':build.sha256(A/'count_selection.json'),
            'safety_selection_script_sha256':build.sha256(ROOT/'select_v8_safe_generator.py'),
            'supported_budget_fidelity':fidelity,
            'safety_guard':'Retain expression <=1.02 baseline and PDS within .05; require neutral false-DE fraction <=.005 in every context/seed and budget-weighted fidelity on reference-call>=10 targets within .02 of baseline.',
            'selection_note':'Calibration-only refinement: weak source targets often have zero or tiny reference DE budgets, making candidate-specific NaN/mean fidelity unstable. A common reference-supported, count-weighted fidelity guard replaces that mean guard. All original raw metric means and the unconstrained count selection are preserved. Final safe winner is frozen before any confirmation scoring.'}
    build.write_json(A/'safe_count_selection.json',record)
    print('V8 safe count winner frozen:',winner,'response',record['response'],'generator',record['generator'],flush=True)

if __name__=='__main__':main()
