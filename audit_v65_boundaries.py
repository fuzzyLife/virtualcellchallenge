"""Audit V65 source exclusion and target-excluded common-template arithmetic."""
import numpy as np,gc
from threadpoolctl import threadpool_limits
import v65_closed_response as m
from v57_pie import sha,write

def main():
 m.check();d=m.old.data();records=[]
 for outer in ('H1','Orion','CD4'):
  b=m.old.bank((outer,));allowed=np.array([m.old.GROUP[str(f)]!=outer for f in m.old.FAMILIES]);assert not b['ok'][~allowed].any()
  before=[tuple(x[allowed].copy() for x in m.source_arrays(b,c)[:2]) for c in (False,True)]
  rows=np.array([m.old.GROUP[str(f)]==outer for f in d['row_families']]);cis=np.unique(d['row_contexts'][rows]);saved={k:d[k][rows].copy() for k in ('truth','errors','valid','knockdown')};profiles=d['profiles'][cis].copy()
  try:
   d['truth'][rows]=17;d['errors'][rows]=.00001;d['valid'][rows]=True;d['knockdown'][rows]=.01;d['profiles'][cis]=1234567;m.old.bank.cache_clear();after=m.old.bank((outer,))
   assert not after['ok'][~allowed].any()
   for c,orig in zip((False,True),before):
    result=m.source_arrays(after,c)
    for p,q in zip(orig,result[:2]):np.testing.assert_array_equal(p,q[allowed])
  finally:
   for k,v in saved.items():d[k][rows]=v
   d['profiles'][cis]=profiles;m.old.bank.cache_clear()
  records.append(dict(held_study=outer,poisoned_rows=int(rows.sum()),allowed_response_and_variance_bitwise_unchanged=True,recipient_support_completely_absent=True));del before,after,b;gc.collect();print('V65 exclusion audit',outer,flush=True)
 write(m.A/'boundary_audit.json',dict(passed=True,checks=records,auditor_sha256=sha(__file__),plan_sha256=sha(m.A/'plan.json'),scope='Recipient treatment, uncertainty, support, knockdown and recipient profile fields poisoned simultaneously; allowed transformed source evidence stays bitwise identical. Separately supplied query controls remain legitimate inference inputs.'))
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
