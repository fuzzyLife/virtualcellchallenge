"""Run one unchanged V88 recipient-prediction job for bounded parallel orchestration."""
import argparse,gc
import numpy as np
from threadpoolctl import threadpool_limits
import v88_mixture_quantiles as m
from v57_pie import read,write,sha
from run_v35_v36 import key

def main(name):
 m.check();dest=m.A/'predictions';dest.mkdir(exist_ok=True);out=dest/(name+'.npz');assert not out.exists()
 if name=='CD4_D3_Rest':
  p=m.ROOT/'experiments/v84_intervention_operator/analysis/donor_predictions.npz';assert sha(p)==read(p.with_suffix('.json'))['sha256']
  with np.load(p,allow_pickle=False) as z:arr={k:z[k].copy() for k in ['query','cp','bulk','v66','genes','targets']}
 else:
  p=m.ROOT/'experiments/v66_assay_strength/analysis/predictions'/(name+'.npz');assert sha(p)==read(p.with_suffix('.json'))['sha256']
  with np.load(p,allow_pickle=False) as z:arr=dict(query=z['query'].copy(),bulk=z['bulk'].copy(),v66=z['assay75'].copy())
 pred,tr=m.predict(arr['query'],arr['v66'],arr['bulk'],key(name));np.savez_compressed(out,**arr,**pred);write(out.with_suffix('.json'),dict(sha256=sha(out),parent_sha256=sha(p),outer=key(name),trace=tr,script_sha256=sha(__file__),plan_sha256=sha(m.A/'plan.json')));print('V88 context complete',name,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('context',choices=['H1','HCT116','HEK293T','CD4_Rest_D1','CD4_Stim48_D2','CD4_D3_Rest']);a=p.parse_args()
 with threadpool_limits(limits=2):main(a.context)
