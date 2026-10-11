"""Apply V15's frozen raw-metric policy without private anchor estimates."""
import argparse,json
import numpy as np,pandas as pd
from scipy.stats import rankdata
from v15_ensemble_methods import A
import build_submission as build
FIELDS=['pds','fidelity','reach','jaccard','nmae']

def compare(frame,name,guards,null):
    base=frame[frame.candidate=='v13'].set_index(['family','seed']);row=frame[frame.candidate==name].set_index(['family','seed']).loc[base.index]
    delta={};complete=True
    for field in FIELDS:
        ok=np.isfinite(base[field])&np.isfinite(row[field]);complete &= bool(ok.all())
        values=(row[field]-base[field])[ok]
        delta[field]=float(values.groupby('family').mean().mean()) if len(values) else None
    den=base.common_expression_denominator;num_base=base.common_expression_numerator;num_row=row.common_expression_numerator
    # The paired numerator ratio equals the expression ratio comparison when the shared denominator resolves.
    relative=(num_row-num_base)/np.maximum(abs(num_base),1e-12)
    expression=float(relative.groupby('family').mean().mean());delta['relative_expression_worsening']=expression
    checks={f:delta[f] is not None and delta[f]>=-guards[f+'_drop'] for f in ['pds','fidelity','reach','jaccard']}
    checks['nmae']=delta['nmae'] is not None and delta['nmae']<=guards['nmae_worsening']
    checks['expression']=np.isfinite(expression) and expression<=guards['relative_expression_worsening']
    for f in ['fidelity','reach']:
        by=(row[f]-base[f]).groupby('family').mean();checks['family_'+f]=bool((by>=-guards['family_'+f+'_drop']).all())
    improvements={f:delta[f] is not None and (delta[f]>0 if f!='nmae' else delta[f]<0) for f in FIELDS};improvements['expression']=expression<0
    checks['three_metrics_improve']=sum(improvements.values())>=3;checks['DE_metric_improves']=any(improvements[f] for f in ['fidelity','reach','jaccard'])
    if 'quantile' in name:
        amount='quantile10' if name.endswith('10') else 'quantile25'
        checks['null_calls']=bool(null.loc[amount,'predicted_calls']<=1.1*null.loc['compact','predicted_calls'])
        if 'shuffled_target25' in frame.candidate.unique():
            shuffled=frame[frame.candidate=='shuffled_target25'].set_index(['family','seed']).loc[base.index]
            wins=[float((row[f]-shuffled[f]).groupby('family').mean().mean())>0 for f in ['fidelity','reach','jaccard']]
            checks['beats_shuffled_distribution']=sum(wins)>=2
    return {'candidate':name,'pass':all(checks.values()),'checks':checks,'changes':delta,'improvements':improvements,'complete_primary_metrics':complete,
            'missing_nmae_panels':[{'family':str(f),'seed':int(s)} for f,s in base.index[~np.isfinite(base.nmae)]],
            'tiny_expression_denominator_panels':int((den<=1e-4).sum()),'scope':'Public partly reused biology; paired absolute expression numerators protect tiny-denominator panels. No statistical or private leaderboard guarantee.'}

def assess(stage):
    plan=json.loads((A/'plan.json').read_text());frame=pd.read_csv(A/f'{stage}_panels.csv');null=pd.read_csv(A/'null_summary.csv').set_index('candidate')
    for p in [stage,'null']:assert json.loads((A/f'{p}_evaluation.json').read_text())['complete']
    if stage=='development':names=['response','v13_quantile25','response_quantile10','response_quantile25']
    else:names=[json.loads((A/'selection.json').read_text())['candidate']]
    comparisons=[compare(frame,n,plan['raw_guards'],null) for n in names]
    passing=[r for r in comparisons if r['pass']]
    rank_scores={}
    if stage=='development' and passing:
        values=np.array([[r['changes'][k] if k!='relative_expression_worsening' else r['changes'][k] for k in FIELDS+['relative_expression_worsening']] for r in passing])
        values[:,4:]*=-1;ranks=np.array([rankdata(-values[:,i],method='average') for i in range(values.shape[1])]).T.mean(1)
        rank_scores={r['candidate']:float(v) for r,v in zip(passing,ranks)};chosen=passing[int(np.argmin(ranks))]['candidate']
    elif stage=='development':chosen='v13'
    else:chosen=names[0]
    result={'stage':stage,'candidate':chosen,'new_candidate_passes':bool(passing),'comparisons':comparisons,'metric_rank_scores':rank_scores,
        'complete_six_metric_verification':bool(passing) and all(r['complete_primary_metrics'] for r in passing),
        'scores_sha256':build.sha256(A/f'{stage}_panels.csv'),'plan_sha256':build.sha256(A/'plan.json'),'null_sha256':build.sha256(A/'null_summary.csv'),
        'scope':'Frozen exploratory public raw policy; reserved algorithm verification on previously inspected biology, not independent biological proof. No retuning.'}
    path=A/('selection.json' if stage=='development' else 'confirmation_result.json')
    if path.exists():raise FileExistsError('Preserve once-only V15 assessment')
    build.write_json(path,result);print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);assess(p.parse_args().stage)
