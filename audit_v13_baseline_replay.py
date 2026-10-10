"""Compare all V13 baseline count summaries against independently generated V12 caches."""
import json
from pathlib import Path
import numpy as np
from v13_methods import ROOT,A
import build_submission as build

def main():
 original=ROOT/'experiments/v12_distribution_calibrated/analysis/harmonized_stage/count_cache'
 checked=[]
 for folder in [A,A/'robust_stage',A/'quality_stage']:
  evaluation=json.loads((folder/'development_evaluation.json').read_text(encoding='utf-8'));assert evaluation['completed']
  for file in sorted((folder/'count_cache').glob('*_development.npz')):
   before=original/file.name;assert before.exists(),before
   with np.load(before) as x,np.load(file) as y:
    for key in ['truth','truth_jack','v11_pred','v11_jack']:np.testing.assert_array_equal(x[key],y[key])
   old=json.loads(before.with_suffix('.json').read_text(encoding='utf-8'));new=json.loads(file.with_suffix('.json').read_text(encoding='utf-8'))
   oldrows={r['reference_mode']:r for r in old['records'] if r['candidate']=='v11'}
   for row in new['records']:
    if row['candidate']!='v11':continue
    for key in ['fidelity','reach','jaccard','nmae','predicted_calls','reference_calls']:
     value=oldrows[row['reference_mode']][key]
     if value is None:self_value=row[key];assert self_value is None
     else:np.testing.assert_allclose(row[key],value,rtol=0,atol=1e-12)
   checked.append({'stage':str(folder.relative_to(A)),'file':file.name})
 build.write_json(A/'baseline_replay_audit.json',{'passed':True,'cases':len(checked),'arrays_per_case':['truth','truth_jack','v11_pred','v11_jack'],'records':checked,'script_sha256':build.sha256(Path(__file__)),'scope':'V11 exact equivalence on reused public splits; not additional independent biological evidence.'})
 print('All baseline profiles, jackknifes, and DE metrics match:',len(checked),flush=True)
if __name__=='__main__':main()
