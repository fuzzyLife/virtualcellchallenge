"""Run the frozen validator with trusted local object-label metadata support only.

The fitter saved its own query genes as NumPy object strings. Numerical arrays,
model, package and frozen validator stay unchanged; only three SHA-pinned query
archives may load their locally generated string objects. All axis assertions run.
"""
from pathlib import Path
import json,runpy
import numpy as np
import build_submission as b
ROOT=b.ROOT;A=ROOT/'experiments/v24_cd4_response/analysis'

def main():
 model=json.loads((A/'reference_model.json').read_text());assert b.sha256(ROOT/'fit_v24_submission.py')==model['dependencies']['fit_v24_submission.py'];assert b.sha256(ROOT/'validate_v24_submission.py')==model['dependencies']['validate_v24_submission.py'];manifest,genes,_=b.read_spec(ROOT);files={str((A/'production_inputs'/(r['context']+'_query.npz')).resolve()):r['query_sha256'] for r in model['context_original_inputs']};original=np.load
 for name,sha in files.items():
  assert b.sha256(Path(name))==sha
  with original(name,allow_pickle=True) as archive:
   labels=archive['genes'];assert labels.shape==(18533,) and all(isinstance(v,(str,np.str_)) for v in labels);np.testing.assert_array_equal(labels.astype(str),genes);assert archive['query_mean'].dtype.kind=='f' and np.isfinite(archive['query_mean']).all()
 def pinned_load(file,*args,**kwargs):
  key=str(Path(file).resolve()) if isinstance(file,(str,Path)) else None
  if key in files:
   assert b.sha256(Path(key))==files[key];kwargs['allow_pickle']=True
  return original(file,*args,**kwargs)
 recovery={'recovery':'Frozen fitter saved locally generated gene labels as object strings; allow NumPy object decoding only for its three checksum-pinned query archives. Every original axis/count/replay assertion remains active. No source/model/count/package values changed.','query_archives':files,'validator_sha256':b.sha256(ROOT/'validate_v24_submission.py'),'fitter_sha256':b.sha256(ROOT/'fit_v24_submission.py'),'wrapper_sha256':b.sha256(Path(__file__))};b.write_json(A/'validation_metadata_recovery.json',recovery)
 try:
  np.load=pinned_load;runpy.run_path(str(ROOT/'validate_v24_submission.py'),run_name='__main__')
 finally:np.load=original
 result=json.loads((A/'independent_validation.json').read_text());assert result['passed'];result['metadata_recovery_sha256']=b.sha256(A/'validation_metadata_recovery.json');b.write_json(A/'independent_validation.json',result)
if __name__=='__main__':main()
