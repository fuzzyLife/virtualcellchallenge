"""Apply V16's prospectively frozen available-metric and null policies."""
import argparse,json
import numpy as np,pandas as pd
from scipy.stats import rankdata
from v16_methods import A,CANDIDATES
from run_v16_experiment import write,load_plan
import build_submission as b
FIELDS=['pds','fidelity','reach','jaccard','nmae']

def compare(frame,name,guards,null):
    base=frame[frame.candidate=='v13'].set_index(['family','seed']);row=frame[frame.candidate==name].set_index(['family','seed']).loc[base.index];delta={};complete=True
    for f in FIELDS:
        ok=np.isfinite(base[f])&np.isfinite(row[f]);complete &= bool(ok.all());values=(row[f]-base[f])[ok];delta[f]=float(values.groupby('family').mean().mean()) if len(values) else None
    relative=(row.common_expression_numerator-base.common_expression_numerator)/np.maximum(abs(base.common_expression_numerator),1e-12);expression=float(relative.groupby('family').mean().mean());delta['relative_expression_worsening']=expression
    checks={f:delta[f] is not None and delta[f]>=-guards[f+'_drop'] for f in FIELDS[:4]};checks['nmae']=delta['nmae'] is not None and delta['nmae']<=guards['nmae_worsening'];checks['expression']=np.isfinite(expression) and expression<=guards['relative_expression_worsening']
    for f in ['fidelity','reach']:checks['family_'+f]=bool(((row[f]-base[f]).groupby('family').mean()>=-guards['family_'+f+'_drop']).all())
    improvements={f:delta[f] is not None and (delta[f]>0 if f!='nmae' else delta[f]<0) for f in FIELDS};improvements['expression']=expression<0;checks['three_metrics_improve']=sum(improvements.values())>=3;checks['DE_metric_improves']=any(improvements[f] for f in ['fidelity','reach','jaccard'])
    kind=CANDIDATES[name][1]
    if kind!='compact':
        checks['null_calls']=bool(null.loc[kind,'predicted_calls']<=1.1*null.loc['compact','predicted_calls']);checks['null_variance_improves']=bool(null.loc[kind,'variance_log_distance']<null.loc['compact','variance_log_distance']);checks['null_detection']=bool(null.loc[kind,'detection_MAE']<=1.1*null.loc['compact','detection_MAE'])
    if name.startswith('calibrated_'):
        uniform=frame[frame.candidate==name.replace('calibrated_','uniform_')].set_index(['family','seed']).loc[base.index]
        checks['beats_amplitude_matched_uniform']=sum(float((row[f]-uniform[f]).groupby('family').mean().mean())>0 for f in ['fidelity','reach','jaccard'])>=2
    return {'candidate':name,'pass':all(checks.values()),'checks':checks,'changes':delta,'improvements':improvements,'complete_primary_metrics':complete,'missing_nmae_panels':[{'family':str(f),'seed':int(s)} for f,s in base.index[~np.isfinite(base.nmae)]],'tiny_expression_denominator_panels':int((base.common_expression_denominator<=1e-4).sum())}

def assess(stage):
    plan=load_plan();file=A/('selection.json' if stage=='development' else 'confirmation_result.json')
    if file.exists():raise FileExistsError('Preserve once-only V16 assessment')
    frame=pd.read_csv(A/f'{stage}_panels.csv');null=pd.read_csv(A/'null_summary.csv').set_index('candidate');assert json.loads((A/f'{stage}_evaluation.json').read_text())['complete']
    selected=json.loads((A/'selection.json').read_text()) if stage!='development' else None
    names=[n for n in CANDIDATES if n!='v13'] if stage=='development' else [selected['candidate']];comparisons=[compare(frame,n,plan['raw_guards'],null) for n in names];passing=[r for r in comparisons if r['pass']];ranks={}
    if stage=='development' and passing:
        values=np.array([[r['changes'][k] for k in FIELDS+['relative_expression_worsening']] for r in passing]);values[:,4:]*=-1;order=np.array([rankdata(-values[:,i]) for i in range(values.shape[1])]).T.mean(1);ranks={r['candidate']:float(v) for r,v in zip(passing,order)};chosen=passing[int(np.argmin(order))]['candidate']
    elif stage=='development':chosen='v13'
    else:chosen=selected['candidate']
    result={'stage':stage,'candidate':chosen,'new_candidate_passes':bool(passing),'comparisons':comparisons,'metric_rank_scores':ranks,'complete_six_metric_verification':bool(passing) and all(r['complete_primary_metrics'] for r in passing),'scope':'Exploratory public partly reused biology. Available raw guards are not a private anchor-scaled score or a statistical rank guarantee.'}
    result.update({'plan_sha256':b.sha256(A/'plan.json'),'scores_sha256':b.sha256(A/f'{stage}_panels.csv'),'null_sha256':b.sha256(A/'null_summary.csv')});write(file,result);print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);assess(p.parse_args().stage)
