"""Independent case/source/cache and baseline-parity audit for V22."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from v22_methods import ROOT,A,load_plan,raw,load_model,write
from fit_v7_response import read_cache
import build_submission as b

def main():
 plan=load_plan();assert b.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')==plan['source_cache_sha256'];data=read_cache();genes=data['genes'];targets=data['targets'];cases=0;baseline_matches=0;cache_count=0
 for ctx in plan['contexts']:
  f=ctx['family'];x,meta,measured,_=raw(f,genes);labels=meta.gene_target.to_numpy(str);assert not set(ctx['left'])&set(ctx['right']);assert (labels[ctx['left']]=='Non-Targeting').all() and (labels[ctx['right']]=='Non-Targeting').all()
  for c in ctx['cells']:
   assert len(set(c['actual']))==len(set(c['template']))==400;assert (labels[c['actual']]==c['target']).all();assert (labels[c['template']]=='Non-Targeting').all();assert set(c['template'])<=set(ctx['left']);assert not set(c['actual'])&(set(ctx['left'])|set(ctx['right']))
   side=A/'development_cache'/f'{f}_{c["target"]}_{c["seed"]}.json';rec=json.loads(side.read_text());assert b.sha256(side.with_suffix('.npz'))==rec['sha256'];cache_count+=1
   base=next(v for v in rec['records'] if v['candidate']=='v16');old=ROOT/'experiments/v17_learned_heterogeneity/analysis/development_cache'/side.name
   prior=json.loads(old.read_text());before=next(v for v in prior['records'] if v['candidate']=='v16');assert base['count_sha256']==before['count_sha256'],'Baseline count mismatch '+str(side);baseline_matches+=1
   assert len(rec['records'])==7 and all(v['actual_cells']==v['prediction_cells']==400 for v in rec['records']);cases+=1
  for c in ctx['oracle']:
   assert len(set(c['calibration']))==len(set(c['actual']))==200;assert not set(c['calibration'])&set(c['actual']);assert (labels[c['actual']]==c['target']).all() and (labels[c['calibration']]==c['target']).all()
  model=load_model(f);record=json.loads((A/'models'/f'{f}.json').read_text());assert f not in record['training_families'];assert set(record['original_inputs'])==set(record['training_families']);assert not np.isin(genes[model['state_axis']],targets).any();assert not np.isin(genes[model['output_axis']],targets).any();np.testing.assert_allclose(model['basis'].T@model['basis'],np.eye(4),atol=1e-6);assert model['interaction'].shape==(300,4,2048);assert np.isfinite(model['interaction']).all() and abs(model['interaction']).max()<=.400001;assert np.isfinite(model['strength_residual']).all();assert (np.diff(model['strength_residual'],axis=1)>=-1e-7).all()
  for family,original in record['original_inputs'].items():
   directory=ROOT/('reference/v5_h1' if family=='H1' else 'reference/v7_sources/'+family);assert original['raw']==json.loads((directory/'extraction.json').read_text());frame=pd.read_csv(directory/'selected_cells.csv');column='target_gene' if family=='H1' else 'gene_target';lab=frame[column].to_numpy(str);assert np.isin(lab[original['control_rows']],['non-targeting','Non-Targeting']).all()
  del x
 oracle=pd.read_csv(A/'oracle_scores.csv');assert len(oracle)==45 and (oracle.actual_cells==200).all() and (oracle.prediction_cells==200).all();null=pd.read_csv(A/'null_scores.csv');assert len(null)==6 and all(v.count_sha256.nunique()==1 for _,v in null.groupby('family'))
 write(A/'independent_experiment_audit.json',{'passed':True,'held_family_cases':cases,'full_count_populations':cases*7,'exact_prior_V16_count_hash_matches':baseline_matches,'cache_hashes_checked':cache_count,'oracle_cases':9,'oracle_populations':45,'oracle_disjoint_calibration_truth_cells':200,'full_axis_null_cases':3,'null_count_identity_verified':True,'source_exclusions_verified':True,'panel_gene_feature_exclusion_verified':True,'orthogonal_state_basis_verified':True,'plan_sha256':b.sha256(A/'plan.json'),'auditor_sha256':b.sha256(Path(__file__)),'limitations':'Public development biology reused; inference factors are statistical proxies and do not identify causal baseline state or guide efficacy. CD4 reserve unopened unless selected.'});print('V22 independent audit passed:',cases,'cases,',baseline_matches,'exact V16 baseline matches',flush=True)
if __name__=='__main__':main()
