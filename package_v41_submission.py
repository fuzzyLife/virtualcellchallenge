"""Lossless bounded-memory reencoding, then unchanged official vcc prep.

Float32 raw integer-valued counts and pre-widened CSR indices avoid SciPy's
large conversion copies above 2**31 nonzeros. Prediction values never change.
"""
from pathlib import Path
import tempfile,json,subprocess,sys,shutil
import h5py,numpy as np
import build_submission as build
import pandas as pd
pd.options.future.infer_string=False
ROOT=build.ROOT
OUT=ROOT/'experiments/v41_identity_contrast'
A=OUT/'analysis'
def load():
 manifest,genes,targets=build.read_spec(ROOT)
 record=json.loads((A/'reference_model.json').read_text())
 assert build.sha256(A/'reference_model_arrays.npz')==record['arrays_sha256']
 assert build.sha256(A/'release_decision.json')==record['selection_sha256']
 for name,h in record['dependencies'].items():assert build.sha256(ROOT/name)==h
 with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as z:
  np.testing.assert_array_equal(z['genes'],genes)
  np.testing.assert_array_equal(z['targets'],targets)
 return manifest,genes,targets,record,None
from vcc.vccfile import validate_vcc,read_vcc_meta

def main():
 manifest,genes,targets,model,effects=load();del effects
 destination=OUT/'prediction.vcc';source=OUT/'prediction.h5ad'
 if destination.exists():raise FileExistsError('Preserve existing upload package')
 verification=json.loads((A/'prediction_validation.json').read_text())
 expected=verification['nnz'];before=source.stat()
 with tempfile.TemporaryDirectory(prefix='v41_official_stage_') as folder:
  scratch=Path(folder).resolve()
  assert scratch.parent==Path(tempfile.gettempdir()).resolve() and scratch.name.startswith('v41_official_stage_')
  staged=scratch/'prediction_float32_int64.h5ad';local_package=scratch/'prediction.vcc'
  counters={}
  with h5py.File(source,'r') as src,h5py.File(staged,'w') as dst:
   for key,value in src.attrs.items():dst.attrs[key]=value
   for key in src:
    if key!='X':src.copy(key,dst)
   group=dst.create_group('X')
   for key,value in src['X'].attrs.items():group.attrs[key]=value
   assert src['X'].attrs['shape'].tolist()==[360000,18533]
   for key,dtype in [('data',np.float32),('indices',np.int64),('indptr',np.int64)]:
    original=src['X/'+key]
    converted=group.create_dataset(key,shape=original.shape,dtype=dtype,chunks=(min(1000000,len(original)),))
    checked=0
    for start in range(0,len(original),4000000):
     block=original[start:start+4000000];value=block.astype(dtype)
     np.testing.assert_array_equal(value,block)
     if key=='data':assert np.issubdtype(block.dtype,np.integer) and (block>0).all() and (block<=1000000).all()
     converted[start:start+len(block)]=value;checked+=len(block)
     if checked%400000000==0:print('Lossless stage',key,f'{checked:,}/{len(original):,}',flush=True)
    counters[key]=checked
  after=source.stat();assert before.st_size==after.st_size and before.st_mtime_ns==after.st_mtime_ns
  assert counters['data']==counters['indices']==expected and counters['indptr']==360001
  source_sha=build.sha256(source)
  build.write_json(A/'memory_safe_stage.json',{'all_values_indices_and_pointers_compared':True,'checked_elements':counters,
   'source_sha256':source_sha,'source_bytes':source.stat().st_size,'staged_bytes':staged.stat().st_size,
   'data_dtype':'float32 with exactly unchanged raw integer values','indices_dtype':'int64','indptr_dtype':'int64',
   'reason':'Uniform lossless staging avoids CSR conversion copies; valid below and above 2**31 nonzeros.',
   'adapter_sha256':build.sha256(Path(__file__))})
  print('Lossless staging verified; invoking unchanged official vcc prep',flush=True)
  result=subprocess.run([sys.executable,'-m','vcc','--json','prep',str(staged),'-g',str(ROOT/'gene_names.csv'),
   '--perts',str(ROOT/'pert_counts.csv'),'--contexts',','.join(manifest['contexts']),
   '--cells-per-pert',str(manifest['cells_per_pert']),'-o',str(local_package)],capture_output=True,text=True,encoding='utf-8',errors='replace')
  (A/'official_prep_memory_stdout.json').write_text(result.stdout,encoding='utf-8')
  (A/'official_prep_memory_stderr.txt').write_text(result.stderr,encoding='utf-8')
  if result.returncode:raise RuntimeError(result.stderr[-2000:])
  validate_vcc(local_package);meta=read_vcc_meta(local_package);assert meta['nnz']==expected
  digest=build.sha256(local_package);pending=destination.with_suffix('.vcc.partial')
  assert pending.resolve().parent==OUT.resolve()
  shutil.copyfile(local_package,pending)
  assert build.sha256(pending)==digest;pending.replace(destination)
  build.write_json(A/'submission_ready.json',{'file':str(destination),'bytes':destination.stat().st_size,'sha256':digest,
   'official_prep_passed':True,'validation_status':model['status'],'fixed_public_development_and_confirmation_guards_passed':model['available_public_guards_passed'],'complete_six_metric_public_verification':model['available_public_guards_passed'],'confirmed_leaderboard_improvement':False,'package_meta':meta,'prediction_validation':verification,'uploaded':False,'leaderboard_score':None,
   'generator_sha256':build.sha256(ROOT/'build_v41_submission.py'),'methods_sha256':build.sha256(ROOT/'v41_identity_contrast.py'),'emitter_sha256':build.sha256(ROOT/'v16_methods.py'),
   'model_sha256':build.sha256(A/'reference_model.json'),'packaging_adapter_sha256':build.sha256(Path(__file__)),
   'memory_safe_stage_sha256':build.sha256(A/'memory_safe_stage.json')})
  print('V41 official package prepared and copied with matching SHA-256',flush=True)
if __name__=='__main__':main()
