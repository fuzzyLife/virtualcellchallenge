"""Apply V14's prospectively declared public raw guards; never retune confirmation."""
import argparse,json
import numpy as np,pandas as pd
from v14_context_methods import A
import build_submission as build
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio']


def assess(stage):
    history=json.loads((A/'metadata_correction.json').read_text())
    plan=json.loads((A/'plan.json').read_text());frame=pd.read_csv(A/f'{stage}_panels.csv')
    summary=pd.read_csv(A/f'{stage}_summary.csv').set_index('candidate')
    candidate=json.loads((A/'screen_result.json').read_text())['advance']
    base=summary.loc['v13'];row=summary.loc[candidate];delta=row[FIELDS]-base[FIELDS]
    primary=frame[frame.reference_mode=='full_independent']
    complete=bool(np.isfinite(primary.loc[primary.candidate.isin(['v13',candidate]),FIELDS]).all().all())
    guards=plan['raw_guards'];family=primary.groupby(['candidate','family'])[FIELDS].mean()
    utility=lambda x:float(np.mean([x.pds,x.fidelity,x.reach,x.jaccard,1-x.nmae,1-x.common_expression_ratio]))
    gains={k:float(delta[k]) for k in FIELDS}
    improvements={k:bool(delta[k]>0 if k in FIELDS[:4] else delta[k]<0) for k in FIELDS}
    checks={'primary_complete':complete,'expression_noninferiority':bool(row.common_expression_ratio<=base.common_expression_ratio*(1+guards['relative_expression_worsening'])),
            'fidelity':bool(delta.fidelity>=-guards['fidelity_drop']),'reach':bool(delta.reach>=-guards['reach_drop']),
            'nmae':bool(delta.nmae<=guards['nmae_worsening']),'pds':bool(delta.pds>=-guards['pds_drop']),
            'at_least_four_metrics_improve':sum(improvements.values())>=4,'DE_yield_or_overlap_improves':any(improvements[k] for k in ['fidelity','reach','jaccard']),
            'utility_improves':utility(row)>utility(base)}
    for field in ['fidelity','reach']:
        checks['family_'+field]=all(family.loc[(candidate,f),field]>=family.loc[('v13',f),field]-guards['family_'+field+'_drop'] for f in family.index.get_level_values('family').unique())
    if stage=='development':checks['beats_shuffled_context_utility']=utility(row)>utility(summary.loc['shuffled2'])
    result={'stage':stage,'candidate':candidate,'promote':all(checks.values()),'checks':checks,'changes':gains,'metric_improvements':improvements,
            'utility':{name:utility(x) for name,x in summary.iterrows()},'development_sha256':build.sha256(A/'development_summary.csv'),
            'panels_sha256':build.sha256(A/f'{stage}_panels.csv'),'plan_sha256':build.sha256(A/'plan.json'),
            'metadata_correction_sha256':build.sha256(A/'metadata_correction.json'),
            'scope':'Prospectively bounded exploratory raw utility and guards on partly reused public biology; reserved verification is not blind biological confirmation. Not private anchored scores or statistical noninferiority. No retuning after verification.'}
    path=A/('selection.json' if stage=='development' else 'confirmation_result.json')
    if path.exists():raise FileExistsError('Preserve frozen selection')
    build.write_json(path,result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);assess(p.parse_args().stage)
