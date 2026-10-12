"""Independent identity, source exclusion, model pooling and metric-file audit."""
import argparse,json
import numpy as np,pandas as pd
import build_submission as b
from v17_methods import ROOT,A,RAW_FAMILIES
from fit_v7_response import raw_inputs,read_cache
from run_v17_experiment import load_plan,load_model,write

def main(stage):
 plan=load_plan();plan_path=A/'plan.json'
 if stage=='confirmation':
  plan_path=A/'confirmation_plan.json';plan=json.loads(plan_path.read_text());assert plan['parent_plan_sha256']==b.sha256(A/'plan.json')
 data=read_cache();checks=[]
 for ctx in plan['contexts']:
  dev={r['target'] for r in ctx['cells'] if r['split']=='development'};verification={r['target'] for r in ctx['cells'] if r['split']=='confirmation'}
  if stage=='confirmation':assert not dev&verification
  family=ctx['family'];left=set(ctx['left']);right=set(ctx['right']);assert not left&right and len(left)>=600
  x,meta,_,_=raw_inputs(family,data['genes']);labels=meta.gene_target.to_numpy(str);assert set(labels[list(left|right)])=={'Non-Targeting'}
  for c in ctx['cells']:
   assert len(c['actual'])==len(c['template'])==400 and len(set(c['actual']))==len(set(c['template']))==400;assert set(c['template'])<=left;assert not set(c['actual'])&(left|right);assert np.all(labels[c['actual']]==c['target'])
  checks.append({'family':family,'control_fit_rows':len(left),'independent_control_rows':len(right),'cases':len(ctx['cells'])})
  del x
 sources={f:pd.read_csv(A/f'{f}_variation_fits.csv') for f in RAW_FAMILIES}
 for f in RAW_FAMILIES+['all']:
  model=load_model(f);training=[v for v in RAW_FAMILIES if v!=f];frame=pd.concat([sources[v] for v in training],ignore_index=True);means=frame.groupby(['family','target_index']).cv2.mean().groupby('target_index').mean();expected=np.zeros(300);expected[means.index.to_numpy(int)]=means.to_numpy();np.testing.assert_array_equal(model['cv2'],expected);np.testing.assert_array_equal(model['generic_cv2'],np.full(300,expected.mean()));np.testing.assert_array_equal(model['shuffled_cv2'],np.random.default_rng(17017).permutation(expected))
 for f,frame in sources.items():
  assert (frame.train_cells+frame.test_cells>=80).all() and (frame.genes>=10).all() and ((frame.cv2>=0)&(frame.cv2<=.5)).all();assert set(frame.fold)=={0,1}
 evaluation=json.loads((A/f'{stage}_evaluation.json').read_text());assert evaluation['complete'] and evaluation['plan_sha256']==b.sha256(plan_path)
 for key,file in [('scores_sha256',f'{stage}_scores.csv'),('panels_sha256',f'{stage}_panels.csv')]:assert evaluation[key]==b.sha256(A/file)
 frame=pd.read_csv(A/f'{stage}_scores.csv');assert not frame.duplicated(['family','target','seed','candidate']).any();assert frame.prediction_cells.eq(400).all()
 for file in (A/(stage+'_cache')).glob('*.json'):
  rec=json.loads(file.read_text());assert rec['sha256']==b.sha256(file.with_suffix('.npz')) and rec['inputs']['plan']==b.sha256(plan_path)
 write(A/f'{stage}_independent_audit.json',{'passed':True,'plan_sha256':b.sha256(plan_path),'auditor_sha256':b.sha256(ROOT/'audit_v17_experiment.py'),'verified_case_identities':checks,'compared_populations':len(frame),'whole_family_exclusions_verified':True,'model_arrays_recomputed_from_source_fits':True,'negative_controls_recomputed_exactly':True,'cache_and_metric_hashes_verified':True,'scope':'Source/exclusion/identity/provenance verification; no independent biological ground truth or hidden anchor score.'});print('V17 independent experiment audit passed',stage,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);main(p.parse_args().stage)
