"""Exact-input memoization of immutable V62 pool outputs; no numerical change."""
import hashlib
import numpy as np
import v62_conserved_response as m

def install():
 if getattr(m.pool,'v62_memoized',False):return
 original=m.pool;cache={}
 def cached(*args,**kwargs):
  h=hashlib.sha256()
  for value in args:
   a=np.ascontiguousarray(value);h.update(str((a.shape,a.dtype.str)).encode());h.update(a.tobytes())
  h.update(repr(sorted(kwargs.items())).encode());key=h.hexdigest()
  if key not in cache:cache[key]=original(*args,**kwargs)
  return cache[key].copy()
 cached.v62_memoized=True;m.pool=cached
