"""Report paired family-level uncertainty and render raw-metric gains for completed stages."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import build_submission as build
ROOT=build.ROOT;A=ROOT/'experiments/v12_distribution_calibrated/analysis'
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio']
def main(folder,stage):
 out=A/folder;meta=json.loads((out/f'{stage}_evaluation.json').read_text());assert meta['completed']
 for suffix in ['scores','panels','summary']:assert build.sha256(out/f'{stage}_{suffix}.csv')==meta[suffix+'_sha256']
 plan=json.loads((out/'validation_plan.json').read_text());assert meta['inputs']['plan']==build.sha256(out/'validation_plan.json')
 sides=list((out/'count_cache').glob(f'*_{stage}.json'));expected=sum(c['split']==stage for ctx in plan['contexts'] for c in ctx['cells']);assert len(sides)==expected
 for side in sides:
  saved=json.loads(side.read_text());assert saved['inputs']==meta['inputs'] and saved['sha256']==build.sha256(side.with_suffix('.npz'))
 panels=pd.read_csv(out/f'{stage}_panels.csv');records=[]
 for mode in ['sampled400','full_independent']:
  family=panels[panels.reference_mode.eq(mode)].groupby(['candidate','family','context'])[FIELDS].mean().groupby(['candidate','family']).mean();families=sorted(set(family.index.get_level_values('family')))
  indices=np.random.default_rng(2026100961).integers(len(families),size=(10000,len(families)))
  for name in sorted(set(panels.candidate)-{'v11'}):
   difference=family.loc[name].loc[families]-family.loc['v11'].loc[families]
   for metric in FIELDS:
    values=difference[metric].to_numpy();samples=values[indices].mean(axis=1);lo,hi=np.quantile(samples,[.025,.975]);records.append({'reference_mode':mode,'candidate':name,'metric':metric,'difference':float(values.mean()),'bootstrap_lower':float(lo),'bootstrap_upper':float(hi),'families':len(families),'worst_family_difference':float(values.max() if metric in FIELDS[-2:] else values.min())})
 result=pd.DataFrame(records);result.to_csv(out/f'{stage}_paired_family_intervals.csv',index=False)
 full=result[result.reference_mode.eq('full_independent')];names=sorted(full.candidate.unique());fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
 for ax,metric in zip(axes.ravel(),FIELDS):
  s=full[full.metric.eq(metric)].set_index('candidate').loc[names];direction=-1 if metric in FIELDS[-2:] else 1;value=direction*s.difference.to_numpy();bounds=np.sort(direction*s[['bootstrap_lower','bootstrap_upper']].to_numpy(),axis=1)
  ax.errorbar(value,np.arange(len(names)),xerr=np.maximum(np.array([value-bounds[:,0],bounds[:,1]-value]),0),fmt='o',capsize=3);ax.axvline(0,color='grey',lw=1);ax.set_yticks(np.arange(len(names)),names);ax.set_title(metric.replace('common_expression_ratio','expression error')+' gain');ax.set_xlabel('Improvement over V11 (higher is better)')
 fig.suptitle('Public held-family count tests: '+folder+' / '+stage+'\nFamily-cluster intervals; few families, conditional on panels and fixed settings',fontsize=12);fig.savefig(out/f'{stage}_metric_gains.png',dpi=160);plt.close(fig)
 build.write_json(out/f'{stage}_evidence_audit.json',{'passed':True,'complete_cache_count':len(sides),'all_cached_profiles_and_input_hashes_verified':True,'evaluation_sha256':build.sha256(out/f'{stage}_evaluation.json'),'script_sha256':build.sha256(Path(__file__)),'interval_scope':'10,000 paired bootstrap resamples of family mean differences; contexts and seeds averaged inside families, not independent replicates. Only 3 development or 2 confirmation families: intervals are unstable, descriptive, conditional on chosen targets/settings, and exclude selection and private-domain uncertainty.'})
 print('Evidence hashes and paired-family intervals complete',folder,stage,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--folder',required=True);p.add_argument('--stage',default='development');a=p.parse_args();main(a.folder,a.stage)
