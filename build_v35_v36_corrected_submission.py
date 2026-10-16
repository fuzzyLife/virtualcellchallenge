"""Plain-string axis serialization repair; scientific predictions unchanged.

Preserves frozen original generator/model files. Never unpickles object arrays:
labels are reread from the canonical competition CSVs, numeric arrays use
allow_pickle=False and are compared exactly after reencoding.
"""
import argparse,json
from pathlib import Path
import numpy as np
import build_v35_v36_submission as b
A=b.A;ROOT=b.ROOT

def repair():
 file=A/'reference_model_arrays.npz';meta=A/'reference_model.json';rec=json.loads(meta.read_text());assert b.b.sha256(file)==rec['arrays_sha256'];initial=A/'reference_model_arrays_object_initial.npz';initial_meta=A/'reference_model_object_initial.json';assert not initial.exists() and not initial_meta.exists()
 _,genes,targets=b.b.read_spec(ROOT)
 with np.load(file,allow_pickle=False) as z:effects=z['effects'].copy();dispersion=z['dispersion'].copy()
 parent=ROOT/'experiments/v24_cd4_response/analysis/reference_model_arrays.npz'
 with np.load(parent,allow_pickle=False) as z:np.testing.assert_array_equal(genes,z['genes']);np.testing.assert_array_equal(targets,z['targets'])
 pending=A/'reference_model_arrays_unicode_pending.npz';np.savez_compressed(pending,effects=effects,dispersion=dispersion,genes=genes.astype(str),targets=targets.astype(str))
 with np.load(pending,allow_pickle=False) as z:
  np.testing.assert_array_equal(z['effects'],effects);np.testing.assert_array_equal(z['dispersion'],dispersion);np.testing.assert_array_equal(z['genes'],genes);np.testing.assert_array_equal(z['targets'],targets);assert z['genes'].dtype.kind=='U' and z['targets'].dtype.kind=='U'
 initial_meta.write_bytes(meta.read_bytes());file.replace(initial);pending.replace(file);old=rec['arrays_sha256'];rec['arrays_sha256']=b.b.sha256(file);rec['axis_serialization_adapter_sha256']=b.b.sha256(Path(__file__));b.b.write_json(meta,rec)
 b.b.write_json(A/'axis_serialization_correction.json',{'initial_metadata_sha256':b.b.sha256(initial_meta),'initial_arrays_sha256':old,'corrected_arrays_sha256':rec['arrays_sha256'],'adapter_sha256':b.b.sha256(Path(__file__)),'numeric_effect_and_dispersion_arrays_exactly_unchanged':True,'gene_and_target_values_exactly_equal_to_canonical_spec_and_V24':True,'object_arrays_unpickled':False,'new_counts_written_before_correction':False,'reason':'Pandas returned object-string axes; frozen strict loader rejected their pickle encoding. Encode canonical labels as NumPy Unicode, retain all numeric predictions and frozen original source. Original arrays/metadata preserved.'});print('Plain-string axes repaired; numeric predictions unchanged',flush=True)

def fit():
 original=np.savez_compressed
 def save(file,*args,**kwargs):
  for key in ['genes','targets']:
   if key in kwargs:kwargs[key]=np.asarray(kwargs[key],str)
  return original(file,*args,**kwargs)
 np.savez_compressed=save
 try:b.fit_production()
 finally:np.savez_compressed=original
 rec=json.loads((A/'reference_model.json').read_text());rec['axis_serialization_adapter_sha256']=b.b.sha256(Path(__file__));b.b.write_json(A/'reference_model.json',rec)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['repair','fit','preflight','generate']);p.add_argument('--workers',type=int,default=4);args=p.parse_args()
 if not 1<=args.workers<=6:p.error('workers1..6')
 if args.stage=='repair':repair()
 elif args.stage=='fit':fit()
 else:b.generate(args.stage=='preflight',args.workers)