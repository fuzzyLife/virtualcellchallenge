"""Pre-transfer correction: absent source genes cannot be zero-valued rank votes."""
import argparse
import numpy as np
import v59_state_composition as m
from run_v59_closed_stage import activate as activate_closure

def activate():
 activate_closure();r=m.read(m.A/'pre_transfer_mask_correction.json')
 assert r['wrapper_sha256']==m.sha(__file__)
 assert r['original_plan_sha256']==m.sha(m.A/'plan.json')
 original=m.mechanism
 def masked(outer,source):
  value=original(outer,source);value['profile']=value['profile'].copy();value['profile'][~value['measured']]=np.nan
  return value
 m.mechanism=masked

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source',choices=m.SOURCES,required=True);a=p.parse_args();activate();m.transfer(a.source)
