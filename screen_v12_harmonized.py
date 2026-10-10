"""Validate Flex source identities and whole-family isolation; screen target-label negatives."""
import json,gc
from pathlib import Path
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
import build_submission as build
from fit_v7_response import read_cache
from v7_response_methods import FAMILIES,summary_score
from v12_harmonized_methods import banks,VARIANTS
ROOT=build.ROOT;A=ROOT/'experiments/v12_distribution_calibrated/analysis/harmonized_stage';D=ROOT/'reference/v12_resources/flex_source'
def main():
 data=read_cache();targets=data['targets'];genes=data['genes'];plan=json.loads((A/'validation_plan.json').read_text())
 reserved={int(np.flatnonzero(targets==c['target'])[0]) for ctx in plan['contexts'] for c in ctx['cells'] if c['split']=='confirmation'}
 meta=pd.read_csv(D/'selected_cells.csv');assert meta.cell_id.is_unique
 ntc=meta.gene_target.eq('Non-Targeting');assert meta.loc[ntc,'guide'].str.match(r'^NO-TARGETg\d+$').all()
 assert set(targets)<=set(meta.gene_target)
 with np.load(D/'summaries.npz') as z:
  assert z['truth'].shape==(300,18533);own=np.array([np.flatnonzero(genes==t)[0] for t in targets])
  assert not z['valid'][np.arange(300),own].any();assert not z['valid'][~z['target_qc_pass']].any()
  assert np.isfinite(z['truth']).all() and np.isfinite(z['errors']).all()
 records=[];checks={}
 for family in FAMILIES:
  bank,training=banks(data,excluded=family);assert family not in training['training_families']
  if family=='K562':
   for name in VARIANTS[1:]+['centered_shuffle']:np.testing.assert_array_equal(bank[name]['consensus20'],bank['v11']['consensus20'])
   checks['whole_K562_exclusion_erases_Flex_exactly']=True
  for ri in np.flatnonzero((data['row_families']==family)&~np.isin(data['row_targets'],list(reserved))):
   ti=data['row_targets'][ri];mask=data['valid'][ri].copy();mask[own[ti]]=False
   for name,b in bank.items():
    metrics=summary_score(b['consensus20'][ti],data['truth'][ri],data['errors'][ri],data['q'][ri],mask)
    records.append({'family':family,'context':int(data['row_contexts'][ri]),'target':targets[ti],'candidate':name,**metrics})
  print('Source checks and negative-control screen complete',family,flush=True);del bank;gc.collect()
 df=pd.DataFrame(records);df.to_csv(A/'source_summary_scores.csv',index=False)
 fields=['mse_ratio','fidelity_proxy','jaccard_proxy','sign_reference_DE','selection_objective']
 summary=df.groupby(['candidate','family','context'])[fields].mean().groupby(['candidate','family']).mean();summary.to_csv(A/'source_summary_family.csv');summary.groupby('candidate').mean().to_csv(A/'source_summary.csv')
 build.write_json(A/'source_checks.json',{'passed':True,'checks':checks,'unique_cells':len(meta),'unique_controls':int(ntc.sum()),'control_guide_identities_verified':True,'all_targets_present':True,'own_gene_excluded':True,'failed_intervention_unavailable':True,'summary_sha256':build.sha256(D/'summaries.npz'),'script_sha256':build.sha256(Path(__file__)),'scope':'Existing public source summaries; normal-q proxies only, not official count tests. Reserved confirmation targets excluded from evaluation.'})
 print(summary.groupby('candidate').mean().to_string(),flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
