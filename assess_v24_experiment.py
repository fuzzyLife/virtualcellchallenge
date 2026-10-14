"""Once-only V24 mean-expert selection and fresh-condition confirmation."""
import argparse,json
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy.stats import rankdata
from v24_methods import A,ELIGIBLE,load_plan,write
import build_submission as b
FIELDS=['pds','fidelity','reach','jaccard','nmae']
def compare(f,name,guards):
 base=f[f.candidate=='v16'].set_index(['family','seed']);row=f[f.candidate==name].set_index(['family','seed']).loc[base.index];changes={};checks={};improved={}
 for k in FIELDS:
  valid=np.isfinite(base[k])&np.isfinite(row[k]);values=(row[k]-base[k])[valid];changes[k]=float(values.groupby('family').mean().mean()) if valid.any() else None;improved[k]=changes[k] is not None and (changes[k]<0 if k=='nmae' else changes[k]>0)
  if k!='nmae':checks[k]=changes[k] is not None and changes[k]>=-guards[k+'_drop']
 changes['relative_expression_worsening']=float(((row.common_expression_numerator-base.common_expression_numerator)/np.maximum(abs(base.common_expression_numerator),1e-12)).groupby('family').mean().mean());improved['expression']=changes['relative_expression_worsening']<0;checks['expression']=changes['relative_expression_worsening']<=guards['relative_expression_worsening'];checks['nmae']=changes['nmae'] is not None and changes['nmae']<=guards['nmae_worsening'];checks['three_metrics_improve']=sum(improved.values())>=3;checks['DE_metric_improves']=any(improved[k] for k in ['fidelity','reach','jaccard'])
 for k in ['fidelity','reach']:checks['family_'+k]=bool(((row[k]-base[k]).groupby('family').mean()>=-guards['family_'+k+'_drop']).all())
 control={'cdf25':'shuffled25','cdf50':'shuffled50','cdf50_de':'shuffled50_de'}[name];other=f[f.candidate==control].set_index(['family','seed']).loc[base.index];checks['beats_alpha_matched_target_shuffle']=sum(float((row[k]-other[k]).groupby('family').mean().mean())>0 for k in ['fidelity','reach','jaccard'])>=2
 return {'candidate':name,'pass':all(checks.values()),'checks':checks,'changes':changes,'improvements':improved,'complete_primary_metrics':bool(np.isfinite(base.nmae).all() and np.isfinite(row.nmae).all()),'missing_nmae_panels':[{'family':str(v[0]),'seed':int(v[1])} for v in base.index[~np.isfinite(base.nmae)]]}
def main(stage):
 plan=load_plan();dest=A/('selection.json' if stage=='development' else 'confirmation_result.json')
 if dest.exists():raise FileExistsError('Preserve once-only assessment')
 evidence=json.loads((A/f'{stage}_evaluation.json').read_text());assert evidence['complete'] and evidence['panels_sha256']==b.sha256(A/f'{stage}_panels.csv');f=pd.read_csv(A/f'{stage}_panels.csv');names=ELIGIBLE if stage=='development' else [json.loads((A/'selection.json').read_text())['candidate']];comp=[compare(f,n,plan['raw_guards']) for n in names];passing=[r for r in comp if r['pass']];chosen='v16';ranks={}
 if passing:
  keys=FIELDS+['relative_expression_worsening'];values=np.array([[r['changes'][k] for k in keys] for r in passing]);values[:,4:]*=-1;rank=np.array([rankdata(-values[:,i]) for i in range(6)]).T.mean(1);ranks={r['candidate']:float(v) for r,v in zip(passing,rank)};chosen=passing[int(np.argmin(rank))]['candidate']
 write(dest,{'candidate':chosen,'new_candidate_passes':bool(passing),'comparisons':comp,'metric_improvement_ranks':ranks,'complete_primary_metrics':bool(passing and next(r for r in passing if r['candidate']==chosen)['complete_primary_metrics']),'plan_sha256':b.sha256(A/'plan.json'),'panels_sha256':b.sha256(A/f'{stage}_panels.csv'),'scope':'Rest raw calibration overlaps public summaries; independent400-cell condition test follows. No private anchor score estimated.'});print(dest.read_text(),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);main(p.parse_args().stage)
