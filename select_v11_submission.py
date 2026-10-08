"""Freeze V11 on development; assess the frozen winner on new count-confirmation targets."""
import argparse,json
import numpy as np,pandas as pd
from v11_methods import ROOT,A,SETTINGS
from select_v8_safe_generator import supported_budget_fidelity
import build_submission as build
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio']


def assess(stage,names):
 scores=pd.read_csv(A/f'{stage}_scores.csv');panels=pd.read_csv(A/f'{stage}_panels.csv')
 summary=pd.read_csv(A/f'{stage}_summary.csv').set_index(['reference_mode','candidate'])
 details={};guards={name:True for name in names};utility={name:0. for name in names}
 for mode,weight in [('full_independent',.7),('sampled400',.3)]:
  frame=scores[scores.reference_mode==mode];full=frame.groupby(['candidate','family','context']).fidelity.mean().groupby(['candidate','family']).mean().groupby('candidate').mean().to_dict()
  keys=['family','context','target','seed'];base=frame[frame.candidate=='v10'][keys+['fidelity']].rename(columns={'fidelity':'baseline'})
  paired=frame.merge(base,on=keys,validate='many_to_one').dropna(subset=['fidelity','baseline'])
  finite=paired.groupby(['candidate','family','context']).fidelity.mean().groupby(['candidate','family']).mean().groupby('candidate').mean().to_dict()
  supported=supported_budget_fidelity(frame)
  pk=['family','context','seed'];allp=panels[(panels.reference_mode==mode)&panels.candidate.isin(names)].copy()
  valid=allp.assign(finite=np.isfinite(allp[FIELDS]).all(axis=1)).groupby(pk).finite.agg(['all','count'])
  use=valid[valid['all']&(valid['count']==len(names))].reset_index()[pk]
  common=allp.merge(use,on=pk,validate='many_to_one')
  common['utility']=(2*(common.pds-.5)+common.fidelity+common.reach+common.jaccard+1-common.nmae.clip(0,3)+1-common.common_expression_ratio.clip(0,4))/6
  util=common.groupby(['candidate','family','context']).utility.mean().groupby(['candidate','family']).mean().groupby('candidate').mean().to_dict()
  families=sorted(common.family.unique().tolist());assert families==['H1','HCT116','HEK293T'],'Missing entire family in complete metrics'
  reference=summary.loc[(mode,'v10')]
  source_fidelity=frame.groupby(['candidate','family','context']).fidelity.mean().groupby(['candidate','family']).mean()
  for name in names:
   row=summary.loc[(mode,name)]
   safe=bool(full[name]>=full['v10']-.01 and finite[name]>=finite['v10']-.01 and
    supported.get(name,0)>=supported.get('v10',0)-.02 and row.pds>=reference.pds-.03 and
    row.common_expression_ratio<=reference.common_expression_ratio*1.10 and
    row.reach>=reference.reach-.03 and row.nmae<=reference.nmae+.015)
   if mode=='full_independent':safe=bool(safe and all(source_fidelity.loc[(name,f)]>=source_fidelity.loc[('v10',f)]-.03 for f in families))
   guards[name]&=safe;utility[name]+=weight*util[name]
  details[mode]={'full_target_mean_fidelity':full,'paired_finite_fidelity':finite,'supported_budget_fidelity':supported,
   'common_panel_utility':util,'complete_panels':len(use),'families_present':families,
   'family_fidelity':{name:{f:float(source_fidelity.loc[(name,f)]) for f in families} for name in names}}
 return summary,guards,utility,details


def main(stage):
 meta=json.loads((A/f'{stage}_evaluation.json').read_text());assert meta['completed']
 for file,key in [(f'{stage}_scores.csv','scores_sha256'),(f'{stage}_panels.csv','panels_sha256'),(f'{stage}_summary.csv','summary_sha256')]:assert build.sha256(A/file)==meta[key]
 assert meta['inputs']['methods']==build.sha256(ROOT/'v11_methods.py') and meta['inputs']['script']==build.sha256(ROOT/'evaluate_v11_submission.py')
 plan=json.loads((A/'validation_plan.json').read_text())
 assert meta['inputs']['plan']==build.sha256(A/'validation_plan.json')
 assert plan['methods_sha256']==build.sha256(ROOT/'v11_methods.py') and plan['script_sha256']==build.sha256(ROOT/'prepare_v11_validation.py')
 keys=['family','context','target','seed','candidate','reference_mode']
 frame=pd.read_csv(A/f'{stage}_scores.csv');assert not frame.duplicated(keys).any()
 expected={(ctx['family'],ctx['context'],cell['target'],cell['seed'],name,mode)
  for ctx in plan['contexts'] for cell in ctx['cells'] if cell['split']==stage
  for name in meta['settings'] for mode in ['sampled400','full_independent']}
 assert set(map(tuple,frame[keys].itertuples(index=False,name=None)))==expected,'Missing or extra evaluation cases'
 expected_cases=sum(cell['split']==stage for ctx in plan['contexts'] for cell in ctx['cells'])
 sides=list((A/'count_cache').glob(f'*_{stage}.json'));assert len(sides)==expected_cases
 for side in sides:
  cached=json.loads(side.read_text());assert cached['inputs']==meta['inputs'] and cached['sha256']==build.sha256(side.with_suffix('.npz'))
 if stage=='development':
  if (A/'count_selection.json').exists():raise FileExistsError('Selection already frozen')
  assert not (A/'confirmation_scores.csv').exists()
  names=[name for name,setting in SETTINGS.items() if setting['eligible']]
  summary,guards,utility,details=assess(stage,names)
  power=json.loads((A/'development_power_validation.json').read_text());assert power['completed']
  assert power['methods_sha256']==build.sha256(ROOT/'v11_methods.py') and power['script_sha256']==build.sha256(ROOT/'audit_v11_cell_power.py')
  assert power['scores_sha256']==build.sha256(A/'development_power_scores.csv') and power['summary_sha256']==build.sha256(A/'development_power_summary.csv')
  for name in names:guards[name]&=power['guards'][name]
  passing=[name for name in names if guards[name]];winner=max(passing,key=lambda name:utility[name])
  build.write_json(A/'count_selection.json',{'candidate':winner,'setting':SETTINGS[winner],'settings':SETTINGS,'confirmation_read':False,
   'passing_candidates':passing,'guards':guards,'weighted_utility':utility,'details':details,
   'development_power_validation_sha256':build.sha256(A/'development_power_validation.json'),
   'development_scores_sha256':build.sha256(A/'development_scores.csv'),'development_summary_sha256':build.sha256(A/'development_summary.csv'),
   'plan_sha256':build.sha256(A/'validation_plan.json'),'methods_sha256':build.sha256(ROOT/'v11_methods.py'),
   'evaluation_script_sha256':build.sha256(ROOT/'evaluate_v11_submission.py'),'selection_script_sha256':build.sha256(ROOT/'select_v11_submission.py'),
   'guard':'Both reference sizes: full/paired FID within .01; supported FID within .02; PDS within .03; expression <=1.10 baseline; reach within .03; NMAE <=baseline+.015. Full-reference per-family FID within .03. Three complete families, identical six-metric panel coverage. V4 compact is diagnostic, not eligible.',
   'utility_policy':'Public raw six-metric proxy, 70% full independent references and 30% sampled400, equal family/context weighting. Not anchored leaderboard score.',
   'neutral_policy':'Bias measured and reported; V10 compact baseline already fails biological control-preservation criterion.'})
  print('V11 development winner',winner,'guards',guards,'utility',utility,flush=True)
 else:
  if (A/'confirmation_result.json').exists():raise FileExistsError('Confirmation already recorded')
  selection=json.loads((A/'count_selection.json').read_text());assert selection['selection_script_sha256']==build.sha256(ROOT/'select_v11_submission.py')
  winner=selection['candidate'];names=list(dict.fromkeys(['v10',winner]));summary,guards,utility,details=assess(stage,names)
  power=json.loads((A/'confirmation_power_validation.json').read_text());assert power['completed']
  assert power['methods_sha256']==build.sha256(ROOT/'v11_methods.py') and power['script_sha256']==build.sha256(ROOT/'audit_v11_cell_power.py')
  assert power['scores_sha256']==build.sha256(A/'confirmation_power_scores.csv') and power['summary_sha256']==build.sha256(A/'confirmation_power_summary.csv')
  guards[winner]&=power['guards'][winner]
  difference=utility[winner]-utility['v10']
  build.write_json(A/'confirmation_result.json',{'candidate':winner,'guard_passed':guards[winner],
   'confirmed_improvement':bool(guards[winner] and difference>0),'weighted_utility_difference':difference,'details':details,
   'power_validation_sha256':build.sha256(A/'confirmation_power_validation.json'),
   'selection_sha256':build.sha256(A/'count_selection.json'),'summary_sha256':build.sha256(A/'confirmation_summary.csv'),
   'policy':'Frozen development winner tested once with V10; no selection or retuning from confirmation. HCT/HEK count targets new, H1 reused biology; all public summaries previously available.'})
  print('V11 confirmation',winner,'guard',guards[winner],'utility difference',difference,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--stage',choices=['development','confirmation'],required=True);main(p.parse_args().stage)
