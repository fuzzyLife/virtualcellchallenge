"""Prospective V20 family-balanced raw guards and negative-control comparisons."""
import argparse,json
import numpy as np,pandas as pd
from scipy.stats import rankdata
from v20_methods import A,load_plan
from run_v20_experiment import write
import build_submission as b
FIELDS=['pds','fidelity','reach','jaccard','nmae']
def compare(frame,name,guards):
 base=frame[frame.candidate=='v16'].set_index(['family','seed']);row=frame[frame.candidate==name].set_index(['family','seed']).loc[base.index];delta={};complete=True
 for f in FIELDS:
  ok=np.isfinite(base[f])&np.isfinite(row[f]);complete &= bool(ok.all());values=(row[f]-base[f])[ok];delta[f]=float(values.groupby('family').mean().mean()) if len(values) else None
 relative=(row.common_expression_numerator-base.common_expression_numerator)/np.maximum(abs(base.common_expression_numerator),1e-12);delta['relative_expression_worsening']=float(relative.groupby('family').mean().mean())
 checks={f:delta[f] is not None and delta[f]>=-guards[f+'_drop'] for f in FIELDS[:4]};checks['nmae']=delta['nmae'] is not None and delta['nmae']<=guards['nmae_worsening'];checks['expression']=np.isfinite(delta['relative_expression_worsening']) and delta['relative_expression_worsening']<=guards['relative_expression_worsening']
 for f in ['fidelity','reach']:checks['family_'+f]=bool(((row[f]-base[f]).groupby('family').mean()>=-guards['family_'+f+'_drop']).all())
 improvements={f:delta[f] is not None and (delta[f]>0 if f!='nmae' else delta[f]<0) for f in FIELDS};improvements['expression']=delta['relative_expression_worsening']<0;checks['three_metrics_improve']=sum(improvements.values())>=3;checks['DE_metric_improves']=any(improvements[f] for f in ['fidelity','reach','jaccard'])
 for control in ['amplitude','shuffled']:
  other=frame[frame.candidate==control].set_index(['family','seed']).loc[base.index];diff={f:float((row[f]-other[f]).groupby('family').mean().mean()) for f in ['fidelity','reach','jaccard']};checks['beats_'+control]=sum(v>0 for v in diff.values())>=2;delta['versus_'+control]=diff
 if name=='stack_joint':
  for control in ['stack_mean','reweight']:
   other=frame[frame.candidate==control].set_index(['family','seed']).loc[base.index];diff={f:float((row[f]-other[f]).groupby('family').mean().mean()) for f in ['fidelity','reach','jaccard']};checks['beats_'+control]=sum(v>0 for v in diff.values())>=2;checks['shape_beats_'+control]=float((row.shape_wasserstein-other.shape_wasserstein).groupby('family').mean().mean())<0;delta['versus_'+control]=diff
 return {'candidate':name,'pass':all(checks.values()),'checks':checks,'changes':delta,'improvements':improvements,'complete_primary_metrics':complete,'missing_nmae_panels':[{'family':str(f),'seed':int(s)} for f,s in base.index[~np.isfinite(base.nmae)]]}
def main(stage):
 plan=load_plan();file=A/('selection.json' if stage=='development' else 'confirmation_result.json')
 if file.exists():raise FileExistsError('Preserve once-only V20 assessment')
 frame=pd.read_csv(A/f'{stage}_panels.csv');ev=json.loads((A/f'{stage}_evaluation.json').read_text());assert ev['complete'] and ev['panels_sha256']==b.sha256(A/f'{stage}_panels.csv')
 audit=json.loads((A/f'{stage}_independent_audit.json').read_text());assert audit['passed'] and audit['plan_sha256']==b.sha256(A/'plan.json')
 names=plan['eligible'] if stage=='development' else [json.loads((A/'selection.json').read_text())['candidate']];comparisons=[compare(frame,n,plan['raw_guards']) for n in names];passing=[v for v in comparisons if v['pass']];ranks={};chosen='v16'
 if passing:
  keys=FIELDS+['relative_expression_worsening'];values=np.array([[r['changes'][k] for k in keys] for r in passing]);values[:,4:]*=-1;order=np.array([rankdata(-values[:,i]) for i in range(6)]).T.mean(1);ranks={r['candidate']:float(v) for r,v in zip(passing,order)};chosen=passing[int(np.argmin(order))]['candidate']
 write(file,{'stage':stage,'candidate':chosen,'new_candidate_passes':bool(passing),'comparisons':comparisons,'metric_rank_scores':ranks,'plan_sha256':b.sha256(A/'plan.json'),'panels_sha256':b.sha256(A/f'{stage}_panels.csv'),'scope':'Public reused biology; family-balanced raw guards; no private anchor-scaled score or guaranteed ranking.'});print(json.dumps(json.loads(file.read_text()),indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);main(p.parse_args().stage)
