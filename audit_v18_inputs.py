"""Recompute V18 inference inputs from frozen NTC identities; never use treated counts."""
from pathlib import Path
import argparse,json
import numpy as np
import build_submission as b
from fit_v7_response import read_cache,raw_inputs
from v18_methods import ROOT,A,load_plan

def main(restore=False):
 plan=load_plan();data=read_cache();train=json.loads((ROOT/'analysis/post_v17_research/PIE_training_meta.json').read_text());axis=np.asarray([g for g in data['genes'] if g in set(train['genes'])],str);pos=np.array([list(data['genes']).index(g) for g in axis]);records=[]
 for ctx in plan['contexts']:
  family=ctx['family'];x,meta,measured,source=raw_inputs(family,data['genes']);ids=np.array(ctx['left']);assert np.all(meta.gene_target.to_numpy(str)[ids]=='Non-Targeting')
  depth=meta.total_counts.to_numpy(float)[ids];y=x[ids].astype(np.float32).multiply((1e4/depth).astype(np.float32)[:,None]).tocsr();y.data=np.log1p(y.data);mean=np.asarray(y.mean(0)).ravel()[pos];file=A/'inputs'/f'{family}.npz'
  if restore and not file.exists():file.parent.mkdir(exist_ok=True);np.savez_compressed(file,genes=axis,competition_genes=data['genes'],targets=data['targets'],ctrl_mean=mean,measured=measured[pos],positions=pos,control_rows=ids)
  assert b.sha256(file)==ctx['pie_controls_sha256']
  with np.load(file,allow_pickle=False) as z:
   for k,v in {'genes':axis,'positions':pos,'control_rows':ids,'ctrl_mean':mean,'targets':data['targets'],'competition_genes':data['genes'],'measured':measured[pos]}.items():np.testing.assert_array_equal(z[k],v)
  for side in (A/'inference').glob(f'{family}_*.json'):
   rec=json.loads(side.read_text());assert rec['inputs_sha256']==b.sha256(file) and rec['plan_sha256']==b.sha256(A/'plan.json')
  records.append({'family':family,'original_NTC_rows':len(ids),'input_sha256':b.sha256(file),'all_arrays_recomputed_exactly':True});print('V18 controls-only input verified',family,flush=True)
 b.write_json(A/'controls_input_independent_audit.json',{'passed':True,'plan_sha256':b.sha256(A/'plan.json'),'records':records,'auditor_sha256':b.sha256(Path(__file__)),'treatment_labels_or_values_used_for_inference':False})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--restore-missing',action='store_true');main(p.parse_args().restore_missing)
