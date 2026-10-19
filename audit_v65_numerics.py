"""Independent scalar reconstruction of a V65 closed-source template coordinate."""
import numpy as np
import v65_closed_response as m
from v57_pie import sha,write

def main():
 rng=np.random.default_rng(202665);f,t,g=3,32,45;own=np.arange(t)
 q=rng.uniform(3,99,(f,g));y=rng.normal(0,.3,(f,t,g));var=rng.uniform(.01,.2,(f,t,g));ok=np.ones((f,t,g),bool)
 # A missing output in one donor must not enter its response-common template.
 ok[0,7,40]=False
 b=dict(profiles=q,y=y,var=var,ok=ok,own=own);actual,v,qp=m.source_arrays(b,True)
 a=np.empty_like(y)
 for fi in range(f):
  q0=q[fi]/sum(q[fi])*1e6
  for ti in range(t):
   factor=2**np.where(ok[fi,ti],y[fi,ti],0);factor[own[ti]]=.1
   q1=q0*factor;q1=q1/sum(q1)*1e6;a[fi,ti]=np.log1p(.05*q1)-np.log1p(.05*q0)
 shared=np.sum(np.sum(ok,axis=2)>0,axis=0)>=2
 checked=0
 for fi in range(f):
  for ti in (0,5,17):
   for gi in (2,18,40):
    ids=[j for j in range(t) if j!=ti and ok[fi,j,gi] and shared[j]]
    template=sum(a[fi,j,gi] for j in ids)/len(ids) if len(ids)>=20 else 0
    np.testing.assert_allclose(actual[fi,ti,gi],a[fi,ti,gi]-template,atol=1e-13);checked+=1
 write(m.A/'numerical_audit.json',dict(passed=True,independent_coordinates_checked=checked,exact_source_own_factor=.1,leave_current_target_out=True,invalid_outputs_excluded=True,auditor_sha256=sha(__file__)))
 print('V65 independent scalar closure/leave-target-out reconstruction passed',checked,'coordinates.')
if __name__=='__main__':main()
