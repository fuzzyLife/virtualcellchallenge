"""Measure V11 dispersion-dependent neutral DE bias on disjoint challenge controls."""
import json,gc
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
from v11_methods import ROOT,A,SEEDS,SETTINGS,emit
from v6_methods import moments,de
import build_submission as build


def main():
 if (A/'decoder_neutral_manifest.json').exists():raise FileExistsError('Neutral diagnostic already frozen')
 manifest,genes,targets=build.read_spec(ROOT);records=[];inputs=[]
 for ci,context in enumerate(manifest['contexts']):
  x,guides,library,_,_,_=build.load_control(ROOT,context,genes,manifest)
  inputs.append({'file':f'context_{context}.h5ad','sha256':build.sha256(ROOT/f'context_{context}.h5ad')})
  for seed in SEEDS:
   rng=np.random.default_rng(seed+ci*113);left=[];right=[]
   for guide in sorted(set(guides)):
    rows=rng.permutation(np.flatnonzero(guides==guide));half=len(rows)//2;left.extend(rows[:half]);right.extend(rows[half:])
   left=np.array(left);right=np.array(right);assert not set(left)&set(right)
   sample=build.guide_stratified_sample(guides[left],library[left],400,rng);template=x[left[sample]]
   _,p,b,_=moments(x[left]);depth=np.asarray(template.sum(axis=1)).ravel();prepared=(depth,p,b)
   reference2000=x[rng.choice(right,2000,replace=False)];reference400=reference2000[:400]
   for name in ['v10','noise35','noise65']:
    count,diag=emit(template,np.arange(400),np.zeros(len(genes)),None,seed+ci,SETTINGS[name],prepared=prepared)
    variance=moments(count)[3].sum()/max(moments(template)[3].sum(),1e-20)
    for label,ref in [('sampled400',reference400),('sampled2000',reference2000)]:
     tested=de(count,ref,genes);calls=int((tested['eligible']&(tested['q']<.05)).sum());n=int(tested['eligible'].sum())
     records.append({'context':context,'seed':seed,'candidate':name,'reference_mode':label,'reference_cells':ref.shape[0],
       'false_de_calls':calls,'tested_genes':n,'false_call_fraction':calls/max(n,1),'variance_ratio':float(variance),**diag})
    print('V11 neutral',context,seed,name,'measured',flush=True)
   pd.DataFrame(records).to_csv(A/'decoder_neutral_scores.csv',index=False)
  del x;gc.collect()
 df=pd.DataFrame(records);summary=df.groupby(['candidate','reference_mode']).agg(max_false_de_calls=('false_de_calls','max'),
  max_false_call_fraction=('false_call_fraction','max'),mean_variance_ratio=('variance_ratio','mean')).reset_index()
 summary.to_csv(A/'decoder_neutral_summary.csv',index=False)
 build.write_json(A/'decoder_neutral_manifest.json',{'diagnostic_only':True,'source_inputs':inputs,
  'scores_sha256':build.sha256(A/'decoder_neutral_scores.csv'),'summary_sha256':build.sha256(A/'decoder_neutral_summary.csv'),
  'script_sha256':build.sha256(ROOT/'benchmark_v11_neutral.py'),'methods_sha256':build.sha256(ROOT/'v11_methods.py'),
  'scope':'Own target omitted, zero response, disjoint provided controls, two seeds, 400 and 2000 references. Synthetic compact populations are biased neutral distributions; not claimed biologically faithful.'})
 print(summary.to_string(index=False),flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
