"""Add lower-cell-count but powered confirmation targets before count confirmation is read."""
import json
from pathlib import Path
import numpy as np,pandas as pd
import build_submission as build
ROOT=build.ROOT;A=ROOT/'experiments/v12_distribution_calibrated/analysis'
def main():
 path=A/'confirmation_power_amendment.json';assert not path.exists()
 plan=json.loads((A/'validation_plan.json').read_text());prior=json.loads((ROOT/'analysis/post_v11_research/future_experiment_plan.json').read_text())
 used=set(prior['excluded_previous_count_targets'])|{c['target'] for ctx in plan['contexts'] for c in ctx['cells'] if c['split']=='development'}
 names=['TAF4','PAPOLA'];assert not set(names)&used
 meta=pd.read_csv(ROOT/'reference/v7_sources/HCT116/selected_cells.csv');meta['sample']=meta['sample'].astype(str)
 targets=pd.read_csv(ROOT/'pert_counts.csv').iloc[:,0].to_numpy(dtype=str);new=[]
 for ci,ctx in enumerate(plan['contexts']):
  if ctx['family']!='HCT116':continue
  for seed in plan['seeds']:
   sample=next(c for c in ctx['cells'] if c['seed']==seed)
   for t in names:
    ids=np.flatnonzero(meta.gene_target.eq(t).to_numpy()&meta['sample'].isin(ctx['batches']).to_numpy());ti=int(np.flatnonzero(targets==t)[0])
    actual=np.random.default_rng(seed+ti+71000).choice(ids,min(400,len(ids)),replace=False).tolist();assert len(actual)>=60
    new.append({'context':ctx['context'],'family':'HCT116','target':t,'seed':seed,'split':'confirmation','actual':actual,'template':sample['template'],'reference':sample['reference']})
 build.write_json(path,{'parent_plan_sha256':build.sha256(A/'validation_plan.json'),'script_sha256':build.sha256(Path(__file__)),'additional_cells':new,
  'reason':'The >=150-cell threshold excluded strong, count-untested HCT116 targets. TAF4 has 89/90 cells, PAPOLA 61/61 by stratum. Include genuine smaller groups to test response recovery; never duplicate cells or claim 400 observed cells.',
  'selection_information':'Previously viewed full-source normal-q power and metadata only, no confirmation count DE/prediction outcomes. TAF4 already reserved in HEK; neither target used in development or earlier count tests.',
  'no_original_confirmation_removed':True,'summary_biology_blind':False})
 print('Confirmation amendment frozen: TAF4 and PAPOLA in both HCT116 strata',flush=True)
if __name__=='__main__':main()
