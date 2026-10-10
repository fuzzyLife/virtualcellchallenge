"""Verify identical V11 populations across independent experiment implementations."""
from pathlib import Path
import json
import numpy as np
import build_submission as build
ROOT=build.ROOT;A=ROOT/'experiments/v12_distribution_calibrated/analysis'
def main():
 reference=A/'flex_stage/count_cache';records=[]
 for folder in ['harmonized_stage','poisson_stage','intervention_stage']:
  shared=sorted(p.name for p in (A/folder/'count_cache').glob('*_development.npz') if (reference/p.name).exists());assert shared
  for name in shared:
   left=reference/name;right=A/folder/'count_cache'/name
   for file in [left,right]:assert json.loads(file.with_suffix('.json').read_text())['sha256']==build.sha256(file)
   with np.load(left) as first,np.load(right) as second:
    for key in ['truth','truth_jack','v11_pred','v11_jack']:np.testing.assert_array_equal(first[key],second[key])
  records.append({'folder':folder,'shared_target_seed_cases':len(shared),'all_truth_and_V11_profiles_jackknives_exact':True})
 build.write_json(A/'cross_stage_baseline_replay.json',{'passed':True,'comparisons':records,'script_sha256':build.sha256(Path(__file__)),'scope':'Exact profile/jackknife equality across shared cells, seeds and baseline implementations; no intervention-count rerun.'})
 print(json.dumps(records,indent=2),flush=True)
if __name__=='__main__':main()
