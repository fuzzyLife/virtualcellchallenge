"""Audit source exclusions and bounded V73 serialized predictions."""
import numpy as np
from threadpoolctl import threadpool_limits
import v73_count_assay as m
from v57_pie import read,write,sha

def main():
 m.check();saved=m.sha;checks=[]
 for excluded in ('H1','Orion','CD4'):
  m.teacher.cache_clear();accessed=[]
  def guarded(path):
   p=m.Path(path)
   if p.parent.name=='sources':
    assert m.study(p.stem)!=excluded,('Recipient evidence touched',p)
    accessed.append(p.stem)
   return saved(path)
  m.sha=guarded
  try:
   for centered in ('raw','assay','center'):
    mean,var,weight=m.teacher((excluded,),centered)
    assert np.isfinite(mean).all() and np.isfinite(var).all() and (var>0).all()
    assert ((weight>=0)&(weight<=1)).all()
   assert set(accessed)=={s for s in m.SOURCES if m.study(s)!=excluded}
  finally:m.sha=saved
  checks.append(dict(excluded_study=excluded,source_access_guard_passed=True))
 d=m.data();own=np.array([list(d['genes']).index(t) for t in d['targets']]);predictions=[]
 for ctx in m.contexts():
  name=ctx['family'];p=m.A/'predictions'/(name+'.npz');assert sha(p)==read(p.with_suffix('.json'))['sha256']
  with np.load(p,allow_pickle=False) as z:
   base=z['v66'];bulk=z['bulk']
   for arm in m.ARMS:
    e=z[arm];assert np.isfinite(e).all() and np.max(abs(e))<=3
    bound=.375 if arm=='assay_center50' else .1875
    assert np.max(abs(e-base))<=bound+1e-6
    assert (e[np.arange(300),own]==0).all()
    np.testing.assert_array_equal(e[:,bulk==0],base[:,bulk==0])
  predictions.append(dict(source=name,prediction_sha256=sha(p)))
 write(m.A/'boundary_audit.json',dict(passed=True,checks=checks,predictions=predictions,all_correction_bounds_verified=True,all_own_effects_zero=True,unsupported_query_genes_unchanged=True,auditor_sha256=sha(__file__)))
 print('V73 source exclusion, correction bounds, own suppression and unsupported-axis audits passed.',flush=True)

if __name__=='__main__':
 with threadpool_limits(limits=2):main()
