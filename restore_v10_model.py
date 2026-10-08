"""Restore the ignored V10 numerical cache from frozen code, inputs and selection.

The tracked manifest is preserved. Restoration must match its exact SHA-256.
"""
import argparse,json,gc
import numpy as np
from threadpoolctl import threadpool_limits
import build_submission as build
from v10_methods import ROOT,A,response_bank
from v8_methods import gated
from fit_v7_response import read_cache

def main(serialization_check=False):
 record=json.loads((A/'reference_model.json').read_text())
 destination=A/'reference_model_arrays.npz'
 pending=A/'reference_model_arrays.restored.npz'
 if pending.exists():raise FileExistsError('Preserve prior reconstruction diagnostic: '+str(pending))
 for file,digest in record['dependencies'].items():assert build.sha256(ROOT/file)==digest,file
 if destination.exists():
  assert build.sha256(destination)==record['model_arrays_sha256']
  if not serialization_check:
   print('V10 cache already matches frozen model');return
  with np.load(destination) as z:arrays={k:z[k].copy() for k in ['genes','targets','effects']}
 else:
  assert build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')==record['source_cache_sha256']
  assert build.sha256(ROOT/'experiments/v4_lineage_response/analysis/source_components.npz')==record['v4_source_components_sha256']
  manifest,genes,targets=build.read_spec(ROOT);data=read_cache()
  np.testing.assert_array_equal(data['genes'],genes);np.testing.assert_array_equal(data['targets'],targets)
  bank,_=response_bank(data);effects=[]
  for input_record in record['context_original_inputs']:assert build.sha256(ROOT/input_record['file'])==input_record['sha256']
  own=np.array([np.flatnonzero(genes==t)[0] for t in targets])
  for context in manifest['contexts']:
   x,_,depth,_,_,_=build.load_control(ROOT,context,genes,manifest)
   query=np.asarray(x.multiply((1e6/depth)[:,None]).mean(axis=0)).ravel()
   values=bank[record['response']] if record['response']!='v4' else response_bank(data,context=context)[0]['v4']
   effects.append(np.array([gated(values[ti].copy(),query,int(own[ti])) for ti in range(300)],dtype=np.float32))
   del x;gc.collect()
  arrays={'genes':np.asarray(genes,dtype=str),'targets':np.asarray(targets,dtype=str),'effects':np.asarray(effects)}
 np.savez_compressed(pending,**arrays)
 assert build.sha256(pending)==record['model_arrays_sha256'],'Exact cache SHA differs; preserve diagnostic and use the recorded runtime/source data'
 if destination.exists():
  assert pending.resolve().parent==A.resolve();pending.unlink()
  print('V10 numerical cache serialization reproduces exact frozen SHA-256')
 else:
  assert destination.resolve().parent==pending.resolve().parent==A.resolve()
  pending.replace(destination);print('V10 numerical cache restored exactly from frozen inputs')
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--serialization-check',action='store_true')
 with threadpool_limits(limits=2):main(p.parse_args().serialization_check)
