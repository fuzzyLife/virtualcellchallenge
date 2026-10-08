"""Freeze V9 calibration winner before reading its confirmation targets."""
import json
import pandas as pd
from v9_methods import ROOT,A,GENERATORS
from select_v8_safe_generator import supported_budget_fidelity
import build_submission as build


def assess(stage):
    summary=pd.read_csv(A/f'count_{stage}_summary.csv').set_index('candidate')
    fidelity=supported_budget_fidelity(pd.read_csv(A/f'count_{stage}_scores.csv'))
    base=summary.loc['legacy'];neutral=json.loads((A/'decoder_neutral_manifest.json').read_text())
    guards={g:bool(g in neutral['passed_generators'] and fidelity.get(g,0)>=fidelity.get('legacy',0)-.02 and
        summary.loc[g,'common_expression_ratio']<=base.common_expression_ratio*1.02 and summary.loc[g,'pds']>=base.pds-.05) for g in GENERATORS}
    return summary,fidelity,guards


def main():
    if (A/'count_selection.json').exists():raise FileExistsError('Frozen calibration exists')
    if (A/'count_confirmation_summary.csv').exists():raise ValueError('Must freeze before confirmation')
    neutral=json.loads((A/'decoder_neutral_manifest.json').read_text())
    assert neutral['methods_sha256']==build.sha256(ROOT/'v9_methods.py')
    assert neutral['scores_sha256']==build.sha256(A/'decoder_neutral_scores.csv')
    for side in (A/'count_cache').glob('*_calibration.json'):
        cached=json.loads(side.read_text())
        assert cached['inputs']['script']==build.sha256(ROOT/'evaluate_v9_submission.py')
        assert cached['inputs']['methods']==build.sha256(ROOT/'v9_methods.py')
        assert cached['inputs']['plan']==build.sha256(A/'validation_plan.json')
        assert cached['sha256']==build.sha256(side.with_suffix('.npz'))
    score,fidelity,guards=assess('calibration');passing=[g for g in GENERATORS if guards[g]]
    winner=max(passing,key=lambda g:score.loc[g,'utility'])
    build.write_json(A/'count_selection.json',{'candidate':winner,'response':'v4_a100','generator':winner,'confirmation_read':False,
        'settings':{g:['v4_a100',g] for g in GENERATORS},'passing_calibration_candidates':passing,
        'supported_budget_fidelity':fidelity,'plan_sha256':build.sha256(A/'validation_plan.json'),
        'calibration_scores_sha256':build.sha256(A/'count_calibration_scores.csv'),
        'calibration_summary_sha256':build.sha256(A/'count_calibration_summary.csv'),
        'neutral_manifest_sha256':build.sha256(A/'decoder_neutral_manifest.json'),
        'script_sha256':build.sha256(ROOT/'select_v9_generator.py'),
        'methods_sha256':build.sha256(ROOT/'v9_methods.py'),'evaluation_script_sha256':build.sha256(ROOT/'evaluate_v9_submission.py'),
        'guard':'Neutral false-DE <=.005 every context/seed, reference-call>=10 budget-weighted fidelity within .02 of legacy, expression <=1.02 legacy, PDS within .05; maximize six-metric public proxy.',
        'utility':'Same equal-family raw six-metric proxy as V8; not official anchored score.'})
    print('V9 calibration winner frozen:',winner,flush=True)

if __name__=='__main__':main()
