"""Freeze V10 from paired calibration evidence; assess once on confirmation."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from v10_methods import ROOT,A,SETTINGS
from select_v8_safe_generator import supported_budget_fidelity
import build_submission as build

def balanced_target_fidelity(frame):
 return frame.groupby(['candidate','family','context']).fidelity.mean().groupby(['candidate','family']).mean().groupby('candidate').mean().to_dict()

def assess(stage):
 scores=pd.read_csv(A/f'count_{stage}_scores.csv')
 panels=pd.read_csv(A/f'count_{stage}_panels.csv')
 summary=pd.read_csv(A/f'count_{stage}_summary.csv').set_index('candidate')
 official_mean=balanced_target_fidelity(scores)
 supported=supported_budget_fidelity(scores)
 keys=['family','context','target','seed']
 ref=scores[scores.candidate=='v4_legacy'][keys+['fidelity']].rename(columns={'fidelity':'base_fidelity'})
 merged=scores.merge(ref,on=keys,validate='many_to_one')
 paired=merged[np.isfinite(merged.fidelity)&np.isfinite(merged.base_fidelity)]
 paired_mean=balanced_target_fidelity(paired)
 common_fields=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio']
 pk=['family','context','seed']
 complete=panels[panels.candidate.isin(SETTINGS)].copy()
 complete['_finite']=np.isfinite(complete[common_fields]).all(axis=1)
 ok=complete.groupby(pk)._finite.agg(['all','count'])
 ok=ok[(ok['all'])&(ok['count']==len(SETTINGS))].reset_index()[pk]
 common=complete.merge(ok,on=pk,validate='many_to_one')
 # Signed whole-gene correction stays visible. No made-up metric for an
 # unresolved expression denominator, and identical panel coverage is used.
 common['paired_utility']=(2*(common.pds-.5)+common.fidelity+common.reach+common.jaccard+
                         1-common.nmae.clip(0,3)+1-common.common_expression_ratio.clip(0,4))/6
 util=common.groupby(['candidate','family','context']).paired_utility.mean().groupby(['candidate','family']).mean().groupby('candidate').mean().to_dict()
 base=summary.loc['v4_legacy'];guards={}
 for name in SETTINGS:
  guards[name]=bool(official_mean.get(name,-1)>=official_mean['v4_legacy']-.01 and
    paired_mean.get(name,-1)>=paired_mean['v4_legacy']-.01 and
    supported.get(name,0)>=supported.get('v4_legacy',0)-.02 and
    summary.loc[name,'pds']>=base.pds-.03 and
    summary.loc[name,'common_expression_ratio']<=base.common_expression_ratio*1.10 and name in util)
 counts=scores.groupby('candidate').agg(finite_fidelity=('fidelity','count'),targets=('target','size'),reference_DE_targets=('reference_calls',lambda v:int((v>0).sum())))
 return summary,guards,{'full_target_mean_fidelity':official_mean,'paired_finite_fidelity':paired_mean,
        'supported_budget_fidelity':supported,'common_panel_utility':util,
        'common_six_metric_panels':len(ok),'missingness':counts.reset_index().to_dict('records')}

def main(stage):
 if stage=='calibration':
  if (A/'count_selection.json').exists():raise FileExistsError('Calibration already frozen')
  if (A/'count_confirmation_scores.csv').exists():raise ValueError('Freeze selection before confirmation')
  plan=json.loads((A/'validation_plan.json').read_text())
  assert plan['script_sha256']==build.sha256(ROOT/'prepare_v10_validation.py')
  assert plan['source_cache_sha256']==build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')
  expected=sum(c['split']=='calibration' for context in plan['contexts'] for c in context['cells'])
  frames=pd.read_csv(A/'count_calibration_scores.csv')
  assert frames[frames.candidate.isin(SETTINGS)].shape[0]==expected*len(SETTINGS)
  inputs={'plan':build.sha256(A/'validation_plan.json'),'script':build.sha256(ROOT/'evaluate_v10_submission.py'),'methods':build.sha256(ROOT/'v10_methods.py')}
  sides=list((A/'count_cache').glob('*_calibration.json'))
  assert len(sides)==expected
  for side in sides:
   cache=json.loads(side.read_text());assert cache['inputs']==inputs
   assert cache['sha256']==build.sha256(side.with_suffix('.npz'))
  neutral=json.loads((A/'decoder_neutral_manifest.json').read_text())
  assert neutral['methods_sha256']==inputs['methods']
  assert neutral['scores_sha256']==build.sha256(A/'decoder_neutral_scores.csv')
  summary,guards,detail=assess(stage)
  passing=[n for n in SETTINGS if guards[n]]
  if not passing:raise ValueError('No candidate meets paired calibration guards')
  winner=max(passing,key=lambda n:detail['common_panel_utility'][n])
  response,generator=SETTINGS[winner]
  build.write_json(A/'count_selection.json',{'candidate':winner,'response':response,'generator':generator,
   'confirmation_read':False,'settings':SETTINGS,'passing_calibration_candidates':passing,'guards':guards,**detail,
   'calibration_scores_sha256':build.sha256(A/'count_calibration_scores.csv'),
   'calibration_summary_sha256':build.sha256(A/'count_calibration_summary.csv'),
   'validation_plan_sha256':build.sha256(A/'validation_plan.json'),
   'methods_sha256':build.sha256(ROOT/'v10_methods.py'),'evaluation_script_sha256':build.sha256(ROOT/'evaluate_v10_submission.py'),
   'selection_script_sha256':build.sha256(Path(__file__)),
   'guard':'Full-target AND paired-finite mean fidelity within .01; supported-budget fidelity within .02; PDS within .03; expression <=1.10 baseline. Equal family/context weighting and common six-metric panel coverage.',
   'neutral_policy':'Neutral DE is measured and reported as distribution-model bias, not suppressed. A synthetic population is a competition prediction, not a biological control reconstruction.'})
  print('V10 frozen winner:',winner,detail,flush=True)
 else:
  if (A/'confirmation_result.json').exists():raise FileExistsError('Confirmation already frozen')
  selection=json.loads((A/'count_selection.json').read_text())
  assert selection['methods_sha256']==build.sha256(ROOT/'v10_methods.py')
  assert selection['evaluation_script_sha256']==build.sha256(ROOT/'evaluate_v10_submission.py')
  assert selection['selection_script_sha256']==build.sha256(Path(__file__))
  summary,guards,detail=assess(stage);winner=selection['candidate']
  difference=detail['common_panel_utility'][winner]-detail['common_panel_utility']['v4_legacy']
  build.write_json(A/'confirmation_result.json',{'candidate':winner,'guard_passed':guards[winner],
   'confirmed_improvement':bool(guards[winner] and difference>0),'utility_difference':difference,**detail,
   'frozen_calibration_selection_sha256':build.sha256(A/'count_selection.json'),
   'summary_sha256':build.sha256(A/'count_confirmation_summary.csv'),
   'policy':'One assessment of the frozen calibration winner. No reselection from confirmation. H1 is reused development; HCT/HEK count targets are new, public source summaries previously available.'})
  print('V10 confirmation:',winner,'difference',difference,'guards',guards[winner],flush=True)
if __name__=='__main__':
 from pathlib import Path
 p=argparse.ArgumentParser();p.add_argument('--stage',choices=['calibration','confirmation'],required=True)
 main(p.parse_args().stage)
