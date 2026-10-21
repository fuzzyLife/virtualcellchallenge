"""Independently check D3 cell splits and poison the excluded response family."""
from collections import Counter
from unittest.mock import patch
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
import v84_intervention_operator as m
import v25_methods as b
import v66_assay_strength as v
from v57_pie import read,write,sha

def main():
 folder=m.ROOT/'experiments/v83_context_expansion/analysis/donor'
 cohort=read(folder/'cohort.json');cells=pd.read_csv(folder/'cells.csv')
 left=set(cohort['left']);right=set(cohort['right'])
 assert len(left)==5006 and len(right)==4994 and not left&right
 assert (cells.iloc[sorted(left|right)].gene_target=='Non-Targeting').all()
 views=[]
 for case in cohort['cases']:
  for view in case['views']:
   actual=view['actual'];template=view['template']
   assert len(actual)==len(set(actual))==len(template)==len(set(template))==400
   assert set(template)<=left and not set(actual)&(left|right)
   assert (cells.iloc[actual].gene_target==case['target']).all()
   assert Counter(cells.iloc[actual].lane_id)==Counter(cells.iloc[template].lane_id)
   views.append(dict(target=case['target'],offset=view['offset'],matched_lanes=True))
 with np.load(m.A/'donor_predictions.npz',allow_pickle=False) as z:
  query=z['query'].copy();bulk=z['bulk'].copy();parent=z['v24'].copy();base=z['v66'].copy();expected=z['program25'].copy()
 def forbidden(*args,**kwargs):raise AssertionError('Excluded CD4 teacher was accessed')
 data=m.data();poison=dict(data);rows=data['row_families']=='CD4';contexts=np.unique(data['row_contexts'][rows])
 assert rows.any()
 for key in ['truth','errors','knockdown','valid']:
  poison[key]=data[key].copy()
  poison[key][rows]=False if key=='valid' else (123.0 if key!='knockdown' else .01)
 poison['profiles']=data['profiles'].copy();poison['profiles'][contexts]=999999.
 with patch.object(b,'make_teacher',forbidden):
  np.testing.assert_array_equal(b.baseline(data,query,('CD4',)),parent)
  np.testing.assert_array_equal(b.baseline(poison,query,('CD4',)),parent)
 # Clear bank cache so neither the original nor poisoned run can bypass recomputation.
 v.old.bank.cache_clear()
 with patch.object(v.old,'data',return_value=poison),patch.object(m,'data',return_value=poison):
  again=v.predict(query,parent,bulk,('CD4',))[0]['assay75']
  np.testing.assert_array_equal(again,base)
  result,_=m.predict(query,again,bulk,('CD4',))
  np.testing.assert_array_equal(result['program25'],expected)
 v.old.bank.cache_clear()
 write(m.A/'donor_exclusion_audit.json',dict(passed=True,controls_fit=len(left),controls_reference=len(right),cohort_views=len(views),views=views,excluded_rows_poisoned=int(rows.sum()),excluded_contexts=contexts.tolist(),CD4_teacher_access_forbidden=True,poisoned_outcomes_uncertainties_doses_masks_and_profiles_leave_V24_V66_V84_exact=True,cohort_sha256=sha(folder/'cohort.json'),prediction_sha256=sha(m.A/'donor_predictions.npz'),auditor_sha256=sha(__file__),scope='Exclusion and metadata integrity only; does not create an independent study or prove hidden-panel generalization.'))
 print('V84 donor exclusion audit passed: 28 lane-matched views, poisoned CD4 evidence has no effect',flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
