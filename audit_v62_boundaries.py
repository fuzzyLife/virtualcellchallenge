"""Independently audit study exclusion and leave-target-out template arithmetic."""
import gc
import numpy as np
from threadpoolctl import threadpool_limits
import v62_conserved_response as m
from v57_pie import sha,write

def main():
 m.check();d=m.data();records=[]
 for outer in ('H1','Orion','CD4'):
  before=m.bank((outer,));rows=np.array([m.GROUP[str(f)]==outer for f in d['row_families']]);cis=np.unique(d['row_contexts'][rows]);saved={k:d[k][rows].copy() for k in ('truth','errors','valid','knockdown')};pp=d['profiles'][cis].copy()
  original={k:v.copy() for k,v in before.items() if k!='profiles'};allowed=np.array([m.GROUP[str(f)]!=outer for f in m.FAMILIES]);profiles=before['profiles'][allowed].copy()
  assert not before['ok'][~allowed].any()
  try:
   d['truth'][rows]=17;d['errors'][rows]=.00001;d['valid'][rows]=True;d['knockdown'][rows]=.01;d['profiles'][cis]=1234567
   m.bank.cache_clear();after=m.bank((outer,))
   for k,v in original.items():np.testing.assert_array_equal(after[k],v,err_msg=outer+':'+k)
   np.testing.assert_array_equal(after['profiles'][allowed],profiles)
  finally:
   for k,v in saved.items():d[k][rows]=v
   d['profiles'][cis]=pp;m.bank.cache_clear()
  records.append(dict(outer=outer,excluded_response_rows=int(rows.sum()),recipient_outcome_poisoning_changes_allowed_predictor_inputs=False));print('V62 boundary audit',outer,flush=True);del before,after,original,saved;gc.collect()
 write(m.A/'boundary_audit.json',dict(passed=True,checks=records,plan_sha256=sha(m.A/'plan.json'),auditor_sha256=sha(__file__),scope='Adversarial recipient outcome/control-profile poisoning leaves all source response ingredients and allowed source profiles bitwise identical. Actual query controls are separately supplied and are legitimate inputs.'))

if __name__=='__main__':
 with threadpool_limits(limits=2):main()
