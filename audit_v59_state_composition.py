"""Independent V59 response-boundary and saved-array audit; never fits a candidate."""
import argparse,platform,importlib.metadata
import numpy as np
import v59_state_composition as m
from run_v59_transfer_stage import activate

def boundary():
 activate();pp=m.plan();assert not any('D3' in s for s in m.SOURCES);results=[]
 for outer in m.OUTERS:
  r=m.read(m.A/'states'/(outer+'.json'));assert all(m.GROUP[s]!=outer for s in r['training_sources'])
  assert not set(r['stable_genes'])&set(m.axes()[1]);assert len(set(r['stable_genes']))==512
 for recipient in m.SOURCES:
  rec=m.read(m.A/'predictions'/(recipient+'.json'));assert all(m.GROUP[s]!=m.GROUP[recipient] for s in rec['training_sources']);assert rec['outer']==m.GROUP[recipient]
  file=m.CACHE/'predictions'/(recipient+'.npz');assert m.sha(file)==rec['sha256'];prior=m.load(file)
  original_mechanism=m.mechanism;original_save=m.save;original_write=m.write
  control_keys={'control_state_counts','control_state_means','profile','measured'};blocked=[];saved=[]
  class Controls(dict):
   def __getitem__(self,key):
    if key not in control_keys:
     blocked.append(key);raise AssertionError('Recipient response field read: '+key)
    return super().__getitem__(key)
  def guarded(outer,source):
   value=original_mechanism(outer,source)
   if source==recipient:
    # Treatment-derived fields become unusable. Only actual control fields survive.
    return Controls({k:value[k] for k in control_keys})
   assert m.GROUP[source]!=m.GROUP[recipient]
   return value
  def compare(path,**arrays):
   assert set(arrays)==set(prior)
   for key in arrays:np.testing.assert_array_equal(arrays[key],prior[key])
   saved.append(True);return rec['sha256']
  try:
   m.mechanism=guarded;m.save=compare;m.write=lambda *args,**kwargs:None
   m.transfer(recipient)
  finally:m.mechanism=original_mechanism;m.save=original_save;m.write=original_write
  assert saved==[True] and not blocked
  own=np.array([np.flatnonzero(prior['genes']==t)[0] for t in prior['targets']])
  for arm in m.ARMS:
   assert np.isfinite(prior[arm]).all();np.testing.assert_array_equal(prior[arm][np.arange(300),own],np.zeros(300))
   if arm in ['joint50','pooled50']:np.testing.assert_array_equal(prior[arm][~prior['available']],prior['v24'][~prior['available']])
  results.append(dict(source=recipient,held_response_fields_unavailable=True,all_prediction_arrays_bitwise_identical=True,own_coordinates_zero=True,unsupported_joint_exact_V24=True))
 m.write(m.A/'boundary_audit.json',dict(passed=True,checks=results,plan_sha256=m.sha(m.A/'plan.json'),correction_sha256=m.sha(m.A/'pre_score_correction.json'),auditor_sha256=m.sha(__file__),D3_opened=False))
 print('V59 whole-study/control-only boundary audit passed',flush=True)

def final():
 pp=m.plan();r=m.read(m.A/'response_assessment.json');count=m.read(m.A/'count_assessment.json');boundary=m.read(m.A/'boundary_audit.json');assert boundary['passed'] and count['complete']
 files=[]
 for path in sorted(m.CACHE.rglob('*.npz')):files.append(dict(file=str(path.relative_to(m.CACHE)),bytes=path.stat().st_size,sha256=m.sha(path)))
 m.write(m.A/'array_manifest.json',dict(files=files,cache=str(m.CACHE),total_bytes=sum(f['bytes'] for f in files),regenerable=True,plan_sha256=m.sha(m.A/'plan.json')))
 runtime={name:importlib.metadata.version(name) for name in ['numpy','scipy','pandas','scikit-learn','anndata','torch']}
 m.write(m.A/'final_audit.json',dict(passed=True,scientific_interpretation='Audit success verifies execution/boundaries, not predictive improvement',release_supported=count['release_supported'],comparisons=count['comparisons'],response_cases=r['cases'],python=platform.python_version(),packages=runtime,response_sha256=m.sha(m.A/'response_cases.csv'),count_sha256=m.sha(m.A/'count_cases.csv'),array_manifest_sha256=m.sha(m.A/'array_manifest.json'),plan_sha256=m.sha(m.A/'plan.json'),correction_sha256=m.sha(m.A/'pre_score_correction.json'),auditor_sha256=m.sha(__file__),D3_opened=False,submission_uploaded=False))
 print('V59 final audit passed; release supported:',count['release_supported'],flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['boundary','final']);a=p.parse_args();{'boundary':boundary,'final':final}[a.stage]()
