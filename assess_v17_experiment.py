"""Apply frozen V17 raw-metric, distribution and negative-control gates."""
import argparse,json
import numpy as np,pandas as pd
from scipy.stats import rankdata
from v17_methods import A
from run_v17_experiment import load_plan,write
import build_submission as b
FIELDS=['pds','fidelity','reach','jaccard','nmae']
def compare(frame,name,guards):
 base=frame[frame.candidate=='v16'].set_index(['family','seed']);row=frame[frame.candidate==name].set_index(['family','seed']).loc[base.index];delta={};complete=True
 for f in FIELDS:
  ok=np.isfinite(base[f])&np.isfinite(row[f]);complete &= bool(ok.all());values=(row[f]-base[f])[ok];delta[f]=float(values.groupby('family').mean().mean()) if len(values) else None
 relative=(row.common_expression_numerator-base.common_expression_numerator)/np.maximum(abs(base.common_expression_numerator),1e-12);delta['relative_expression_worsening']=float(relative.groupby('family').mean().mean());delta['shape_wasserstein']=float((row.shape_wasserstein-base.shape_wasserstein).groupby('family').mean().mean())
 checks={f:delta[f] is not None and delta[f]>=-guards[f+'_drop'] for f in FIELDS[:4]};checks['nmae']=delta['nmae'] is not None and delta['nmae']<=guards['nmae_worsening'];checks['expression']=np.isfinite(delta['relative_expression_worsening']) and delta['relative_expression_worsening']<=guards['relative_expression_worsening']
 for f in ['fidelity','reach']:checks['family_'+f]=bool(((row[f]-base[f]).groupby('family').mean()>=-guards['family_'+f+'_drop']).all())
 improvements={f:delta[f] is not None and (delta[f]>0 if f!='nmae' else delta[f]<0) for f in FIELDS};improvements['expression']=delta['relative_expression_worsening']<0;checks['three_metrics_improve']=sum(improvements.values())>=3;checks['DE_metric_improves']=any(improvements[f] for f in ['fidelity','reach','jaccard']);checks['response_distribution_improves']=delta['shape_wasserstein']<0
 for control in ['generic','shuffled']:
  other=frame[frame.candidate==control].set_index(['family','seed']).loc[base.index];checks['beats_'+control]=sum(float((row[f]-other[f]).groupby('family').mean().mean())>0 for f in ['fidelity','reach','jaccard'])>=2
 return {'candidate':name,'pass':all(checks.values()),'checks':checks,'changes':delta,'improvements':improvements,'complete_primary_metrics':complete,'missing_nmae_panels':[{'family':str(f),'seed':int(s)} for f,s in base.index[~np.isfinite(base.nmae)]]}
def main(stage):
 plan=load_plan();file=A/('selection.json' if stage=='development' else 'confirmation_result.json')
 if file.exists():raise FileExistsError('Preserve once-only V17 assessment')
 frame=pd.read_csv(A/f'{stage}_panels.csv');r=json.loads((A/f'{stage}_evaluation.json').read_text());assert r['complete'] and r['panels_sha256']==b.sha256(A/f'{stage}_panels.csv')
 null=pd.read_csv(A/'null_scores.csv');assert null.null_bitwise_v16_identity.all()
 names=plan['eligible_finalists'] if stage=='development' else [json.loads((A/'selection.json').read_text())['candidate']];comparisons=[compare(frame,n,plan['raw_guards']) for n in names];passing=[v for v in comparisons if v['pass']];ranks={}
 chosen='v16'
 if passing:
  keys=FIELDS+['relative_expression_worsening'];values=np.array([[r['changes'][k] for k in keys] for r in passing]);values[:,4:]*=-1;order=np.array([rankdata(-values[:,i]) for i in range(6)]).T.mean(1);ranks={r['candidate']:float(v) for r,v in zip(passing,order)};chosen=passing[int(np.argmin(order))]['candidate']
 write(file,{'stage':stage,'candidate':chosen,'new_candidate_passes':bool(passing),'comparisons':comparisons,'metric_rank_scores':ranks,'plan_sha256':b.sha256(A/'plan.json'),'panels_sha256':b.sha256(A/f'{stage}_panels.csv'),'null_sha256':b.sha256(A/'null_scores.csv'),'scope':'Public reused biology; family-balanced raw guards and joint-response shape, not private anchor-scaled score or guaranteed ranking.'});print(json.dumps(json.loads(file.read_text()),indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);main(p.parse_args().stage)
