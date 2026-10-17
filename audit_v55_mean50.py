"""Independent V55 arithmetic, production and final delivery audit."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v55_mean50_diagnostic/analysis'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
 return h.hexdigest()
def main(stage):
 if stage=='production':
  rec=json.loads((A/'reference_model.json').read_text());assert sha(A/'reference_model_arrays.npz')==rec['arrays_sha256']
  with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as x,np.load(ROOT/'experiments/v24_cd4_response/analysis/reference_model_arrays.npz',allow_pickle=False) as b,np.load(ROOT/'experiments/v36_population_views/analysis/reference_model_arrays.npz',allow_pickle=False) as r:
   for k in ('genes','targets'):np.testing.assert_array_equal(x[k],b[k]);np.testing.assert_array_equal(x[k],r[k])
   expected=.5*b['effects'].astype(float)+.5*r['effects'].astype(float)
   np.testing.assert_allclose(x['effects'],expected,atol=2.4e-7,rtol=0)
   fallback=b['effects']==r['effects'];np.testing.assert_array_equal(x['effects'][fallback],b['effects'][fallback])
   assert np.count_nonzero(x['effects']-b['effects'])>0
  out={'passed':True,'independent_float64_blend_replayed':True,'fallback_exact':True,'arrays_sha256':rec['arrays_sha256'],'auditor_sha256':sha(Path(__file__))}
  (A/'production_effect_audit.json').write_text(json.dumps(out,indent=2)+'\n')
 elif stage=='delivery':
  import audit_v35_v36_delivery as old
  old.A=A;old.OUT=A.parent;old.main()
 print('V55 audit passed',stage,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['production','delivery']);a=p.parse_args();main(a.stage)
