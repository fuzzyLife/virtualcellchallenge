"""Repair V10 validation adapter depth indexing without changing frozen production.

Preserves original validation records. Reassesses the frozen winner only.
"""
import argparse,json,gc
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
from v10_methods import ROOT,A,SEEDS,response_bank,emit
from v8_methods import gated,panel_metrics
from v6_methods import moments,de,direction_metrics
from fit_v7_response import read_cache,raw_inputs
from audit_transfer_validation import profile_jackknife_chunked
from evaluate_atlas_transfer import load_atlas_model
from select_v8_safe_generator import supported_budget_fidelity
import build_submission as build

def production():
 provenance=pd.read_csv(A/'prediction_control_rows.csv')
 np.testing.assert_array_equal(provenance.original_total_counts,provenance.predicted_total_counts)
 assert len(provenance)==360000
 independent=json.loads((A/'independent_validation.json').read_text())
 assert independent['passed'] and independent['verified_rows']==360000
 build.write_json(A/'production_depth_audit.json',{'passed':True,'verified_cells':len(provenance),
  'every_generated_cell_preserves_its_recorded_sampled_control_depth':True,
  'provenance_sha256':build.sha256(A/'prediction_control_rows.csv'),
  'package_sha256':independent['package_sha256'],'script_sha256':build.sha256(ROOT/'audit_v10_validation_depths.py')})
 print('All 360,000 production depths match their sampled controls',flush=True)

def main():
 if (A/'corrected_depth_validation.json').exists():raise FileExistsError('Corrected diagnostic already recorded')
 plan=json.loads((A/'validation_plan.json').read_text());selection=json.loads((A/'count_selection.json').read_text())
 data=read_cache();genes=data['genes'];targets=data['targets'];atlas,_=load_atlas_model()
 winner=selection['candidate'];settings={'v4_legacy':('v4','legacy'),winner:(selection['response'],selection['generator'])}
 records=[];panels=[];depth_mismatches=[]
 for family in ['H1','HCT116','HEK293T']:
  bank,_=response_bank(data,excluded=family)
  x,meta,measured,source=raw_inputs(family,genes)
  for context in [c for c in plan['contexts'] if c['family']==family]:
   population=x[context['left']];depth,p,b,_=moments(population)
   prepared=(depth,p,b,np.asarray(population.mean(axis=0)).ravel());query=p*1e6
   for seed in SEEDS:
    cells=[c for c in context['cells'] if c['split']=='confirmation' and c['seed']==seed]
    templates=x[cells[0]['template']];template_depth=np.asarray(templates.sum(axis=1)).ravel()
    correct_prepared=(template_depth,p,b,prepared[3])
    depth_mismatches.append({'context':context['context'],'seed':seed,'mismatched_old_validation_depths':int((depth[:400]!=template_depth).sum()),
     'old_validation_mean_depth':float(depth[:400].mean()),'intended_template_mean_depth':float(template_depth.mean())})
    refs={'sampled400':x[cells[0]['reference']],'full_independent':x[context['right']]}
    reference_profiles={};reference_jacks={}
    for label,ref in refs.items():reference_profiles[label],reference_jacks[label]=profile_jackknife_chunked(ref)
    actual_profiles=[];actual_jacks=[];pred_profiles={n:[] for n in settings};pred_jacks={n:[] for n in settings}
    for cell in cells:
     target=cell['target'];ti=int(np.flatnonzero(targets==target)[0]);own=int(np.flatnonzero(genes==target)[0])
     actual=x[cell['actual']];tp,tj=profile_jackknife_chunked(actual);actual_profiles.append(tp);actual_jacks.append(tj)
     truths={label:de(actual,ref,genes) for label,ref in refs.items()}
     for name,(response,generator) in settings.items():
      effect=gated(bank[response][ti].copy(),query,own)
      # Legacy uses full-population mean depth for its existing zero-activation
      # rule. Compact/selective use the actual template depths for row constraints.
      use_prepared=prepared if generator=='legacy' else correct_prepared
      count,diag=emit(templates,np.arange(400),effect,own,seed+ti,generator,prepared=use_prepared,atlas=atlas)
      if generator!='legacy':np.testing.assert_array_equal(np.asarray(count.sum(axis=1)).ravel(),template_depth)
      pp,pj=profile_jackknife_chunked(count);pred_profiles[name].append(pp);pred_jacks[name].append(pj)
      for label,ref in refs.items():
       metrics=direction_metrics(de(count,ref,genes),truths[label],genes,own)
       records.append({'family':family,'context':context['context'],'target':target,'seed':seed,
        'candidate':name,'reference_mode':label,'reference_cells':ref.shape[0],'actual_cells':actual.shape[0],**metrics})
      del count
     print('Corrected depths checked',context['context'],seed,target,flush=True)
    for label in refs:
     for name in settings:
      score=panel_metrics(np.array(pred_profiles[name]),np.array(actual_profiles),reference_profiles[label],np.array(pred_jacks[name]),
       np.array(actual_jacks),reference_jacks[label],genes,[c['target'] for c in cells],targets)
      subset=pd.DataFrame([r for r in records if r['candidate']==name and r['context']==context['context'] and r['seed']==seed and r['reference_mode']==label])
      panels.append({'family':family,'context':context['context'],'seed':seed,'candidate':name,'reference_mode':label,
       **subset.mean(numeric_only=True).to_dict(),**score})
    pd.DataFrame(records).to_csv(A/'corrected_depth_scores.csv',index=False)
    pd.DataFrame(panels).to_csv(A/'corrected_depth_panels.csv',index=False)
    pd.DataFrame(depth_mismatches).to_csv(A/'validation_depth_mismatches.csv',index=False)
  del x;gc.collect()
 frame=pd.DataFrame(records);pf=pd.DataFrame(panels)
 fields=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio']
 summary=pf.groupby(['reference_mode','candidate','family','context'])[fields].mean().groupby(['reference_mode','candidate','family']).mean().groupby(['reference_mode','candidate']).mean()
 summary.reset_index().to_csv(A/'corrected_depth_summary.csv',index=False)
 details={}
 for mode in refs:
  ff=frame[frame.reference_mode==mode];base=summary.loc[(mode,'v4_legacy')];row=summary.loc[(mode,winner)]
  full=ff.groupby(['candidate','family','context']).fidelity.mean().groupby(['candidate','family']).mean().groupby('candidate').mean().to_dict()
  keys=['family','context','target','seed']
  baseline=ff[ff.candidate=='v4_legacy'][keys+['fidelity']].rename(columns={'fidelity':'baseline'})
  matched=ff.merge(baseline,on=keys).dropna(subset=['fidelity','baseline'])
  paired=matched.groupby(['candidate','family','context']).fidelity.mean().groupby(['candidate','family']).mean().groupby('candidate').mean().to_dict()
  supported=supported_budget_fidelity(ff)
  complete=pf[pf.reference_mode==mode].copy();pk=['family','context','seed']
  valid=complete.assign(finite=np.isfinite(complete[fields]).all(axis=1)).groupby(pk).finite.all()
  use=complete.merge(valid[valid].reset_index()[pk],on=pk)
  use['utility']=(2*(use.pds-.5)+use.fidelity+use.reach+use.jaccard+1-use.nmae.clip(0,3)+1-use.common_expression_ratio.clip(0,4))/6
  utility=use.groupby(['candidate','family','context']).utility.mean().groupby(['candidate','family']).mean().groupby('candidate').mean().to_dict()
  passed=bool(full[winner]>=full['v4_legacy']-.01 and paired[winner]>=paired['v4_legacy']-.01 and
   supported.get(winner,0)>=supported.get('v4_legacy',0)-.02 and row.pds>=base.pds-.03 and
   row.common_expression_ratio<=base.common_expression_ratio*1.10)
  details[mode]={'full_target_mean_fidelity':full,'paired_finite_fidelity':paired,'supported_budget_fidelity':supported,
   'paired_six_metric_utility':utility,'guard_passed':passed,'positive_proxy_difference':bool(utility[winner]>utility['v4_legacy'])}
 build.write_json(A/'corrected_depth_validation.json',{'candidate':winner,'settings_unchanged':True,'count_depths_exactly_match_actual_templates':True,
  'diagnostics':details,'selection_sha256':build.sha256(A/'count_selection.json'),'plan_sha256':build.sha256(A/'validation_plan.json'),
  'script_sha256':build.sha256(ROOT/'audit_v10_validation_depths.py'),'methods_sha256':build.sha256(ROOT/'v10_methods.py'),
  'scores_sha256':build.sha256(A/'corrected_depth_scores.csv'),'summary_sha256':build.sha256(A/'corrected_depth_summary.csv'),
  'scope':'Corrected validation adapter only. Original records preserved, no parameter tuning or reselection. Previously confirmed public targets are development, not new blind evidence. Both seeds and 400/full independent reference controls.'})
 print(summary.to_string(),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--production-only',action='store_true');args=p.parse_args()
 with threadpool_limits(limits=2):
  if args.production_only:production()
  else:main()
