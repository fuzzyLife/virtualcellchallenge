"""Verify pre-generation truth and count evaluation used identical real populations."""
import numpy as np
from v58_pie import *

def main():
 p=plan();rec=read(A/'raw_head_evaluation.json');assert rec['complete'] and rec['initial_run_before_count_generation']
 checked=[]
 for f in sorted((A/'head_truth').glob('*.json')):
  h=read(f);c=read(A/'counts'/f.name)
  for key in ['actual_rows','reference_rows','template_rows']:assert h[key]==c[key]
  before=CACHE/'head_truth'/f.with_suffix('.npz').name;after=CACHE/'count_statistics'/f.with_suffix('.npz').name
  assert sha(before)==h['sha256'] and sha(after)==c['statistics_sha256']
  with np.load(before,allow_pickle=False) as a,np.load(after,allow_pickle=False) as b:
   np.testing.assert_array_equal(a['q'],b['truth_q']);np.testing.assert_array_equal(a['lfc'],b['truth_lfc'])
  checked.append(f.stem)
 assert len(checked)==44
 write(A/'truth_population_replay_audit.json',{'passed':True,'cases':len(checked),'case_names':checked,'plan_sha256':sha(A/'plan.json'),'statement':'All pre-generation DE/LFC checks and count tests used identical distinct real cells, references and templates, with bitwise identical truth q-values and LFCs.'})
 print('Truth/population replay passed:',len(checked),'cases',flush=True)

if __name__=='__main__':main()
