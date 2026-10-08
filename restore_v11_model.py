"""Restore an ignored V11 model cache and require its exact frozen SHA-256."""
import json,gc
import numpy as np
from v11_methods import ROOT,A,response_bank,response_for
from fit_v7_response import read_cache
from v6_methods import moments
import build_submission as build


def main():
 path=A/'reference_model_arrays.npz'
 if path.exists():raise FileExistsError('Preserve the existing cache')
 model=json.loads((A/'reference_model.json').read_text())
 for name,digest in model['dependencies'].items():assert build.sha256(ROOT/name)==digest,name
 data=read_cache();assert build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')==model['source_cache_sha256']
 manifest,genes,targets=build.read_spec(ROOT);bank,_=response_bank(data);effects=[]
 for record in model['context_original_inputs']:
  assert build.sha256(ROOT/record['file'])==record['sha256']
  x,_,_,_,_,_=build.load_control(ROOT,record['context'],genes,manifest);query=moments(x)[1]*1e6
  effects.append(np.array([response_for(bank,model['setting'],query,int(np.flatnonzero(genes==t)[0]),ti) for ti,t in enumerate(targets)],dtype=np.float32))
  del x;gc.collect()
 pending=A/'restored_model.partial.npz'
 np.savez_compressed(pending,genes=np.asarray(genes,dtype=str),targets=np.asarray(targets,dtype=str),effects=np.asarray(effects))
 assert build.sha256(pending)==model['model_arrays_sha256'],'Exact model restoration failed; preserve frozen manifest and investigate runtime versions'
 pending.replace(path);print('V11 cache restored with exact frozen SHA-256',flush=True)
if __name__=='__main__':main()
