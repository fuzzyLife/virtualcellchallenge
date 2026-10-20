"""Audit V70 source exclusion and serialized teacher evidence."""
import numpy as np
from threadpoolctl import threadpool_limits
import v70_count_teacher as m
from v57_pie import read,write,sha

def main():
 m.check();d=m.data();records=[]
 for name in m.SOURCES:
  path=m.A/'sources'/(name+'.npz');r=read(path.with_suffix('.json'));assert sha(path)==r['sha256']
  assert r['plan_sha256']==sha(m.A/'plan.json')
  with np.load(path,allow_pickle=False) as z:
   y=z['effect'];v=z['variance'];ok=z['valid'];np.testing.assert_array_equal(z['genes'],d['genes']);np.testing.assert_array_equal(z['targets'],d['targets'])
   assert y.shape==v.shape==ok.shape==(300,18533)
   assert np.isfinite(y).all() and np.isfinite(v).all() and (v>0).all() and np.max(abs(y))<=3
   assert not ok[np.arange(300),[list(d['genes']).index(t) for t in d['targets']]].any()
  records.append(dict(source=name,sha256=r['sha256'],targets=len(r['targets']),supported_coordinates=int(ok.sum())))
 # Fail on any attempt to read excluded evidence, before hashing or array loading.
 saved=m.sha
 for excluded in ('H1','Orion','CD4'):
  touched=[]
  def checked(path):
   p=m.Path(path)
   if p.parent.name=='sources':
    assert m.study(p.stem)!=excluded,('Recipient treatment source accessed',p)
    touched.append(p.stem)
   return saved(path)
  m.sha=checked
  try:
   y,w=m.pooled((excluded,));sy,sw=m.pooled((excluded,),True)
   np.testing.assert_array_equal(sy,np.roll(y,1,axis=0));np.testing.assert_array_equal(sw,np.roll(w,1,axis=0))
   assert np.isfinite(y).all() and ((w>=0)&(w<=1)).all()
   assert set(touched)=={s for s in m.SOURCES if m.study(s)!=excluded}
  finally:m.sha=saved
 write(m.A/'boundary_audit.json',dict(passed=True,sources=records,excluded_studies=['H1','Orion','CD4'],excluded_source_access_blocked=True,cyclic_shuffle_exact=True,recipient_outcomes_in_prediction_interface=False,script_sha256=sha(__file__)))
 print('V70 source validity and whole-study exclusion audit passed.',flush=True)

if __name__=='__main__':
 with threadpool_limits(limits=2):main()
