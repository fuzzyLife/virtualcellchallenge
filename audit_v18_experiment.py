"""Independent identity, source exclusion and per-case arithmetic checks for V18."""
import argparse,json
import numpy as np,pandas as pd
import build_submission as b
from v18_methods import ROOT,A,load_plan,load_predictions,blend
from fit_v7_response import raw_inputs,read_cache
from run_v18_experiment import write

def main(stage):
 plan=load_plan();data=read_cache();checks=[]
 for ctx in plan['contexts']:
  family=ctx['family'];left=set(ctx['left']);right=set(ctx['right']);assert not left&right
  dev={r['target'] for r in ctx['cells'] if r['split']=='development'};ver={r['target'] for r in ctx['cells'] if r['split']=='confirmation'};assert not dev&ver
  x,meta,measured,_=raw_inputs(family,data['genes']);labels=meta.gene_target.to_numpy(str);assert set(labels[list(left|right)])=={'Non-Targeting'}
  for c in ctx['cells']:
   assert len(c['actual'])==len(c['template'])==400 and len(set(c['actual']))==len(set(c['template']))==400;assert set(c['template'])<=left;assert not set(c['actual'])&(left|right);assert np.all(labels[c['actual']]==c['target'])
  expert=load_predictions(family,stage);assert np.array_equal(expert['genes'],data['genes'][expert['positions']]);assert np.array_equal(expert['measured'],measured[expert['positions']])
  bank=json.loads((ROOT/'experiments/v16_control_preserving/analysis/source_banks_cache'/f'{family}.json').read_text());assert family not in bank['training_families']
  checks.append({'family':family,'control_fit_rows':len(left),'independent_control_rows':len(right),'development_target_count':len(dev),'reserved_target_count':len(ver)});del x
 ev=json.loads((A/f'{stage}_evaluation.json').read_text());assert ev['complete'] and ev['plan_sha256']==b.sha256(A/'plan.json')
 for key,file in [('scores_sha256',f'{stage}_scores.csv'),('panels_sha256',f'{stage}_panels.csv')]:assert ev[key]==b.sha256(A/file)
 frame=pd.read_csv(A/f'{stage}_scores.csv');assert not frame.duplicated(['family','target','seed','candidate']).any();assert frame.prediction_cells.eq(400).all()
 for file in (A/(stage+'_cache')).glob('*.json'):
  rec=json.loads(file.read_text());assert rec['sha256']==b.sha256(file.with_suffix('.npz')) and rec['inputs']['plan']==b.sha256(A/'plan.json')
 # Exact paired baseline population digests must match the already-tested V17 V16 comparator.
 if stage=='development':
  old=pd.read_csv(ROOT/'experiments/v17_learned_heterogeneity/analysis/development_scores.csv');old=old[old.candidate=='v16'].set_index(['family','target','seed']);current=frame[frame.candidate=='v16'].set_index(['family','target','seed']);assert (old.loc[current.index].count_sha256==current.count_sha256).all()
 write(A/f'{stage}_independent_audit.json',{'passed':True,'plan_sha256':b.sha256(A/'plan.json'),'auditor_sha256':b.sha256(ROOT/'audit_v18_experiment.py'),'identities':checks,'compared_populations':len(frame),'whole_family_exclusion_verified':True,'V16_count_digests_match_previous_development':stage=='development','cache_metric_hashes_verified':True});print('V18 independent audit passed',stage,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);main(p.parse_args().stage)
