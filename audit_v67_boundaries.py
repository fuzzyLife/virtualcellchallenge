"""Poison all recipient-study evidence and require unchanged V67 source models."""
import gc
import numpy as np
from threadpoolctl import threadpool_limits
import v67_intervention_network as m
from v57_pie import write,sha


def main():
    m.check();d=m.old.data();records=[]
    for outer in ('H1','Orion','CD4'):
        excluded=tuple([outer]);original=m.bank(excluded)
        snapshots=[{k:v.copy() if isinstance(v,np.ndarray) else v for k,v in s.items()} for s in original]
        families=[str(f) for f in m.old.FAMILIES if m.old.GROUP[str(f)]==outer]
        rows=np.isin(d['row_families'],families);cis=np.unique(d['row_contexts'][rows])
        saved={k:d[k][rows].copy() for k in ('truth','errors','valid','knockdown')}
        pp=d['profiles'][cis].copy()
        d['truth'][rows]=2.7;d['errors'][rows]=.001;d['valid'][rows]=True;d['knockdown'][rows]=.01;d['profiles'][cis]=1e7
        m.bank.cache_clear();m.old.bank.cache_clear()
        poisoned=m.bank(excluded)
        assert len(snapshots)==len(poisoned)
        for before,after in zip(snapshots,poisoned):
            assert m.old.GROUP[after['family']]!=outer
            for field in ('w','axis','own','profile','present','responses','variance'):
                np.testing.assert_array_equal(before[field],after[field])
        for k,v in saved.items():d[k][rows]=v
        d['profiles'][cis]=pp
        m.bank.cache_clear();m.old.bank.cache_clear();del original,snapshots,poisoned;gc.collect()
        records.append(dict(outer=outer,recipient_rows_poisoned=int(rows.sum()),passed=True))
        print('V67 recipient poisoning passed',outer,flush=True)
    write(m.A/'boundary_audit.json',dict(passed=True,records=records,plan_sha256=sha(m.A/'plan.json'),auditor_sha256=sha(__file__)))


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
