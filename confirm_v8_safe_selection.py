"""Assess the pre-confirmation frozen safe V8 candidate, without reselection."""
import json
import pandas as pd
import build_submission as build
from v8_methods import A
from select_v8_safe_generator import supported_budget_fidelity


def main():
    selection=json.loads((A/'safe_count_selection.json').read_text())
    assert selection['original_count_selection_sha256']==build.sha256(A/'count_selection.json')
    original=json.loads((A/'confirmation_result.json').read_text())
    assert original['summary_sha256']==build.sha256(A/'count_confirmation_summary.csv')
    table=pd.read_csv(A/'count_confirmation_summary.csv').set_index('candidate')
    candidate=selection['candidate'];base=table.loc['baseline_legacy'];row=table.loc[candidate]
    fidelity=supported_budget_fidelity(pd.read_csv(A/'count_confirmation_scores.csv'))
    guard=bool(fidelity.get(candidate,0)>=fidelity.get('baseline_legacy',0)-.02 and row.common_expression_ratio<=1.02*base.common_expression_ratio and row.pds>=base.pds-.05)
    passed=bool(guard and row.utility>base.utility)
    build.write_json(A/'safe_confirmation_result.json',{'candidate':candidate,'confirmed_improvement':passed,'guard_passed':guard,
        'supported_budget_fidelity':fidelity,'utility_difference':float(row.utility-base.utility),'summary_sha256':build.sha256(A/'count_confirmation_summary.csv'),
        'frozen_calibration_selection_sha256':build.sha256(A/'safe_count_selection.json'),
        'raw_confirmation_sha256':build.sha256(A/'confirmation_result.json'),
        'policy':'Final safe calibration winner frozen before any confirmation scoring; confirmation assessed once with no reselection. Hidden competition performance remains unknown.'})
    print('V8 safe confirmation:',candidate,'passed',passed,'utility difference',float(row.utility-base.utility),flush=True)

if __name__=='__main__':main()
