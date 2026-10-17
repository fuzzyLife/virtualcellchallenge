"""Execution-only per-thread cache adapter for the frozen V46 emitter."""
import threading
class PerThreadCache:
 def __init__(self):self.local=threading.local()
 def values(self):
  if not hasattr(self.local,'data'):self.local.data={}
  return self.local.data
 def get(self,*args):return self.values().get(*args)
 def __getitem__(self,key):return self.values()[key]
 def clear(self):self.values().clear()
 def update(self,**kwargs):self.values().update(**kwargs)

