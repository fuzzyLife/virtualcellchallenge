"""Replay the published response pipeline under the current runtime before final adaptation."""
from pathlib import Path
import argparse,gc
import numpy as np
from threadpoolctl import threadpool_limits
import build_submission as b
from v57_pie import read,write,sha
R=Path(__file__).resolve().parent;A=R/'experiments/final_v91/analysis'
def main(context):
 from fit_v7_response import read_cache
 import v24_methods as v24
 import v66_assay_strength as v66
 import v90_diagnostic_release as v90
 d=read_cache();teacher=v24.make_teacher(True);ci=['A','B','C'].index(context)
 with np.load(R/'experiments/v66_assay_strength/analysis/production_inputs'/(context+'_query.npz'),allow_pickle=False) as z:q=z['query'].copy();bulk=z['bulk'].copy()
 records=[]
 for version,path in [(24,'v24_cd4_response'),(66,'v66_assay_strength'),(90,'v90_interval_diagnostic')]:
  with np.load(R/'experiments'/path/'analysis/reference_model_arrays.npz',allow_pickle=False) as z:expected=z['effects'][ci].copy()
  if version==24:actual=v24.predict(d,teacher,q,production=True)[0]['cdf50_de']
  elif version==66:actual=v66.predict(q,base,bulk)[0]['assay75']
  else:actual=v90.predict(q,base,bulk,'production',controls=False)[0]['interval25']
  error=float(np.max(abs(actual-expected)));np.testing.assert_allclose(actual,expected,rtol=0,atol=1e-5);records.append(dict(version=version,maximum_absolute_LFC_error=error,bitwise_equal=bool(np.array_equal(actual,expected)),tolerance=1e-5));base=expected;print('Validation response replay',context,version,error,flush=True);gc.collect()
 write(A/('validation_replay_'+context+'.json'),dict(passed=True,context=context,comparisons=records,script_sha256=sha(Path(__file__)),scope='Numerical response replay under the working final runtime; no new biological accuracy evidence.'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('context',choices=['A','B','C']);args=p.parse_args()
 with threadpool_limits(limits=2):main(args.context)
