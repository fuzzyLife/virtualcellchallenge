"""Reconstruct all production effects through the ORIGINAL transfer implementation.

The production kernel is a separate implementation. This adapter supplies only
query control fields and V24 baseline arrays, captures outputs without writing,
and compares every one of the 900 effect vectors bitwise.
"""
import numpy as np
import v61_pooled_state_response as v

def verify():
 v.activate();m=v.m;record=m.read(m.A/'reference_model.json');file=m.A/'reference_model_arrays.npz';assert m.sha(file)==record['model_arrays_sha256'];arrays=m.load(file);genes,targets=m.axes();np.testing.assert_array_equal(arrays['genes'],genes);np.testing.assert_array_equal(arrays['targets'],targets);checks=[]
 for ci,context in enumerate(['A','B','C']):
  query=m.A/'production_inputs'/(context+'_query.npz');q=m.load(query);q={k:q[k] for k in ['control_state_counts','control_state_means','profile','measured']};original_mechanism=m.mechanism;original_load=np.load;original_save=m.save;original_write=m.write;saved=[]
  def mechanism(outer,source):return q if source==context else original_mechanism(outer,source)
  def redirect(path,*args,**kwargs):
   name=str(path).replace('\\','/')
   if name.endswith('/inputs/'+context+'.npz') and 'v58' in name:return original_load(query,*args,**kwargs)
   return original_load(path,*args,**kwargs)
  def compare(path,**values):
   np.testing.assert_array_equal(values['joint50'],arrays['effects'][ci]);np.testing.assert_array_equal(values['v24'],arrays['parent_effects'][ci]);saved.append(True);return 'verification-only'
  try:
   m.mechanism=mechanism;np.load=redirect;m.save=compare;m.write=lambda *args,**kwargs:None;m.transfer(context,outer='all')
  finally:m.mechanism=original_mechanism;np.load=original_load;m.save=original_save;m.write=original_write
  assert saved==[True];checks.append(dict(context=context,coordinates=int(arrays['effects'][ci].size),bitwise_equal=True));print('V61 all production effects independently replayed',context,flush=True)
 m.write(m.A/'production_effect_replay.json',dict(passed=True,checks=checks,independent_from_production_kernel=True,original_transfer_code_sha256=m.sha(m.ROOT/'v59_state_composition.py'),verifier_sha256=m.sha(__file__),model_sha256=m.sha(m.A/'reference_model.json'),D3_opened=False))
if __name__=='__main__':verify()
