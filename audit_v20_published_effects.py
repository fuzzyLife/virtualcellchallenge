"""Audit response magnitude changes after publication, without hidden treatment labels."""
from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v20_stack_transfer/analysis'
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 new=A/'reference_model_arrays.npz';old=ROOT/'experiments/v16_control_preserving/analysis/reference_model_arrays.npz'
 with np.load(new,allow_pickle=False) as z:current=z['effects'].copy();genes=z['genes'].copy();targets=z['targets'].copy()
 with np.load(old,allow_pickle=False) as z:baseline=z['effects'].copy();np.testing.assert_array_equal(z['genes'],genes);np.testing.assert_array_equal(z['targets'],targets)
 rows=[];inputs={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in [new,old]}
 for ci,context in enumerate(['A','B','C']):
  p=A/'production_inputs'/(context+'_query.npz');inputs[str(p.relative_to(ROOT)).replace('\\','/')]=sha(p)
  with np.load(p,allow_pickle=False) as z:q=z['query_mean'].copy()
  eligible=np.broadcast_to(q>5,(300,len(genes))).copy()
  for ti,target in enumerate(targets):eligible[ti,np.flatnonzero(genes==target)[0]]=False
  before=np.abs(baseline[ci]).astype(float);after=np.abs(current[ci]).astype(float);diff=after-before;changed=eligible&(current[ci]!=baseline[ci]);row={'context':context,'eligible_changed_coordinates':int(changed.sum()),'magnitude_increased_coordinates':int((changed&(diff>0)).sum()),'magnitude_reduced_coordinates':int((changed&(diff<0)).sum()),'opposite_sign_coordinates':int((changed&(current[ci]*baseline[ci]<0)).sum()),'previously_zero_activated_coordinates':int((changed&(baseline[ci]==0)).sum()),'eligible_absolute_LFC_sum_v16':float(before[eligible].sum()),'eligible_absolute_LFC_sum_v20':float(after[eligible].sum()),'relative_LFC_mass_change':float(after[eligible].sum()/before[eligible].sum()-1)};assert row['magnitude_increased_coordinates']+row['magnitude_reduced_coordinates']==row['eligible_changed_coordinates'];rows.append(row)
 output={'scope':'Pre-emission LFC magnitude on query-control CPM>5, own targets excluded. Includes non-DE genes; no private treated truth, significance or causal attribution.','contexts':rows,'input_sha256':inputs,'auditor_sha256':sha(Path(__file__))};(A/'published_effect_magnitude_diagnostic.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8');print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
