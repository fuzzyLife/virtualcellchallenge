"""Paired fidelity uncertainty, source tradeoffs, and exported scientific figure."""
import json
import numpy as np,pandas as pd
from v10_methods import ROOT,A
import build_submission as build

def main():
 selected=json.loads((A/'count_selection.json').read_text())['candidate']
 reports=[];family_rows=[]
 if (ROOT/'experiments/v9_activation_preserving/analysis/reference_model_arrays.npz').exists():
  with np.load(A/'reference_model_arrays.npz') as n,np.load(ROOT/'experiments/v9_activation_preserving/analysis/reference_model_arrays.npz') as v:
   np.testing.assert_array_equal(n['genes'],v['genes']);np.testing.assert_array_equal(n['targets'],v['targets'])
   new=n['effects'];old=v['effects'];changes=[]
   for i,c in enumerate('ABC'):
    changes.append({'context':c,'changed_effects':int((new[i]!=old[i]).sum()),'total_effects':int(new[i].size),
     'new_nonzero_effects':int((new[i]!=0).sum()),'new_mean_absolute_effect':float(abs(new[i]).mean()),
     'old_mean_absolute_effect':float(abs(old[i]).mean()),'mean_absolute_change':float(abs(new[i]-old[i]).mean())})
  build.write_json(A/'response_change.json',{'new_arrays_sha256':build.sha256(A/'reference_model_arrays.npz'),
   'v4_effect_arrays_sha256':build.sha256(ROOT/'experiments/v9_activation_preserving/analysis/reference_model_arrays.npz'),
   'contexts':changes,'scope':'New response estimates compared with verified unchanged V4/V7/V9 effects.'})
 else:
  build.write_json(A/'response_change.json',{'available':False,'reason':'Optional historical V9 ignored cache is absent; it is not required for V10 generation or model integrity.'})
 for stage in ['calibration','confirmation']:
  if stage=='confirmation' and (A/'corrected_depth_scores.csv').exists():
   d=pd.read_csv(A/'corrected_depth_scores.csv');d=d[d.reference_mode=='sampled400']
  else:d=pd.read_csv(A/f'count_{stage}_scores.csv')
  keys=['family','context','target','seed']
  base=d[d.candidate=='v4_legacy'][keys+['fidelity']].rename(columns={'fidelity':'baseline'})
  chosen=d[d.candidate==selected][keys+['fidelity']]
  pairs=chosen.merge(base,on=keys,validate='one_to_one').dropna(subset=['fidelity','baseline'])
  pairs['difference']=pairs.fidelity-pairs.baseline
  # Seeds and technical strata are repeated measurements, not new biology.
  targets=pairs.groupby(['family','target'])[['fidelity','baseline','difference']].mean().reset_index()
  rng=np.random.default_rng(20271009);boot=[]
  families=targets.family.unique()
  for _ in range(2000):
   sample=[]
   for family in rng.choice(families,len(families),replace=True):
    values=targets[targets.family==family].difference.to_numpy()
    sample.append(float(rng.choice(values,len(values),replace=True).mean()))
   boot.append(np.mean(sample))
  ci=np.quantile(boot,[.025,.975])
  reports.append({'stage':stage,'candidate':selected,'paired_target_family_groups':len(targets),
   'hierarchical_bootstrap_fidelity_difference_CI95':ci.tolist(),
   'balanced_target_fidelity_difference':float(targets.groupby('family').difference.mean().mean()),
   'bootstrap_replicates':2000,'scope':'Public paired cases, clustered by biological family and target. Three source families cannot estimate unseen-lineage uncertainty reliably; interval is descriptive, not a hidden-score confidence interval.'})
  for family,g in targets.groupby('family'):
   family_rows.append({'stage':stage,'family':family,'candidate':selected,'unique_targets':len(g),
     'baseline_fidelity':float(g.baseline.mean()),'candidate_fidelity':float(g.fidelity.mean()),'difference':float(g.difference.mean())})
 pd.DataFrame(family_rows).to_csv(A/'paired_family_fidelity.csv',index=False)
 build.write_json(A/'paired_uncertainty.json',{'analyses':reports,'source_score_hashes':{stage:build.sha256(A/('corrected_depth_scores.csv' if stage=='confirmation' and (A/'corrected_depth_scores.csv').exists() else f'count_{stage}_scores.csv')) for stage in ['calibration','confirmation']},'confirmation_scope':'Corrected-depth sampled400 diagnostic; calibration remains the original historical snapshot','script_sha256':build.sha256(ROOT/'analyze_v10_tradeoffs.py')})
 if (A/'corrected_depth_scores.csv').exists():
  corrected=pd.read_csv(A/'corrected_depth_scores.csv');family_checks=[]
  for (mode,family),g in corrected.groupby(['reference_mode','family']):
   values=g.groupby(['candidate','context']).fidelity.mean().groupby('candidate').mean()
   family_checks.append({'reference_mode':mode,'family':family,'unique_targets':g.target.nunique(),
    'baseline_fidelity':float(values['v4_legacy']),'candidate_fidelity':float(values[selected]),
    'difference':float(values[selected]-values['v4_legacy'])})
  pd.DataFrame(family_checks).to_csv(A/'corrected_depth_family_fidelity.csv',index=False)
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axes=plt.subplots(1,3,figsize=(12,3.6))
 metrics=['fidelity','common_expression_ratio','pds'];labels=['Direction fidelity (higher better)','Expression ratio (lower better)','Discrimination (higher better)']
 for ax,metric,label in zip(axes,metrics,labels):
  for i,stage in enumerate(['calibration','confirmation']):
   if (A/'corrected_depth_summary.csv').exists():
    frame=pd.read_csv(A/'corrected_depth_summary.csv');mode='sampled400' if i==0 else 'full_independent';frame=frame[frame.reference_mode==mode].set_index('candidate')
   else:frame=pd.read_csv(A/f'count_{stage}_summary.csv').set_index('candidate')
   ax.bar(i-.16,frame.loc['v4_legacy',metric],width=.3,color='#687782',label='V4 reference' if i==0 else None)
   ax.bar(i+.16,frame.loc[selected,metric],width=.3,color='#176b99',label='V10' if i==0 else None)
  ax.set_xticks([0,1],['400 controls','Full independent']);ax.set_title(label,fontsize=10)
  ax.spines[['top','right']].set_visible(False)
 axes[0].legend(frameon=False,fontsize=8)
 fig.suptitle('Corrected source validation: raw proxies, not leaderboard scores',fontsize=11)
 fig.tight_layout();fig.savefig(A/'validation_comparison.png',dpi=180);plt.close(fig)
 print(pd.DataFrame(family_rows).to_string(index=False),flush=True)
if __name__=='__main__':main()
