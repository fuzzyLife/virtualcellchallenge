"""Pre-score V59 CPM-closure correction; frozen original numerical code preserved."""
import argparse
from functools import lru_cache
import numpy as np
from threadpoolctl import threadpool_limits
import v59_state_composition as m
ORIGINAL_RECONSTRUCT=m.reconstruct

def closed_reconstruct(nc,cm,tilt,within=None):
 value=ORIGINAL_RECONSTRUCT(nc,cm,tilt,within)
 assert np.isfinite(value).all() and (value>=0).all() and value.sum()>0
 return value*(1e6/value.sum())

def activate():
 r=m.read(m.A/'pre_score_correction.json');assert r['original_plan_sha256']==m.sha(m.A/'plan.json')
 assert r['wrapper_sha256']==m.sha(__file__)
 # Cache immutable source matrices within this process, no arithmetic changes.
 m.source_data=lru_cache(maxsize=5)(m.source_data)
 m.reconstruct=closed_reconstruct

def check():
 cm=np.array([[1.,9.],[9.,1.],[4.,6.],[6.,4.]])*1e5;n=np.array([40.,40.,10.,10.]);zero=np.zeros(4)
 np.testing.assert_allclose(closed_reconstruct(n,cm,zero),(n/n.sum())@cm,rtol=1e-14)
 value=closed_reconstruct(n,cm,np.array([.2,-.2,.1,-.1]),np.tile([.5,-.5],(4,1)))
 np.testing.assert_allclose(value.sum(),1e6,rtol=1e-14);assert (value>0).all()
 print('V59 CPM closure and zero-response checks passed',flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['check','states','mechanisms','transfer','evaluate']);p.add_argument('--source',choices=m.SOURCES);p.add_argument('--outer',choices=(*m.OUTERS,'all'));a=p.parse_args();activate()
 with threadpool_limits(limits=2):
  if a.stage=='check':check()
  elif a.stage=='states':m.fit_states(a.outer)
  elif a.stage=='mechanisms':m.mechanisms(a.outer,a.source)
  elif a.stage=='transfer':m.transfer(a.source)
  else:m.evaluate()
