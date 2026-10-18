"""Run the frozen V57 auditor with immutable JSON reads memoized.

This changes I/O only. Original auditor, formulas, samples and plan stay byte-identical.
"""
from functools import lru_cache
import argparse
import audit_v57_pie as frozen
from v57_pie import A,ROOT,sha,write

def main(stage):
 frozen.read=lru_cache(maxsize=None)(frozen.read)
 {'evidence':frozen.evidence,'artifacts':frozen.artifacts}[stage]()
 write(A/'auditor_execution_note.json',{'original_auditor_sha256':sha(ROOT/'audit_v57_pie.py'),'wrapper_sha256':sha(ROOT/'audit_v57_pie_cached.py'),'change':'Memoize reads of immutable JSON inputs; no formula, sample, reference or experiment-plan change','plan_preserved':True})

if __name__=='__main__':
 q=argparse.ArgumentParser();q.add_argument('stage',choices=['evidence','artifacts']);main(q.parse_args().stage)
